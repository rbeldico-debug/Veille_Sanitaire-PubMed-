import gspread
from typing import Set, List
from src.config import settings, PROJECT_DIR
from src.domain.models import EnrichedArticle, SearchConfig
from src.logger import logger
import re


class GoogleSheetsRepository:
    """
    Gère la persistance et la lecture de configuration dans Google Sheets.
    """
    CONFIG_TAB_NAME = "_ADMIN_CONFIG"

    def __init__(self):
        creds_path = PROJECT_DIR / settings.GOOGLE_CREDENTIALS_FILE
        if not creds_path.exists():
            raise FileNotFoundError(f"Fichier credentials introuvable : {creds_path}")

        self.gc = gspread.service_account(filename=str(creds_path))

        try:
            self.sh = self.gc.open(settings.GOOGLE_SHEET_NAME)
            # Par défaut, on n'a pas d'onglet actif, il sera défini par l'orchestrateur
            self.worksheet = None
        except gspread.SpreadsheetNotFound:
            raise ValueError(f"Google Sheet '{settings.GOOGLE_SHEET_NAME}' introuvable.")

    def get_configs(self) -> List[SearchConfig]:
        max_allowed_results = 20  # Sécurité : On ne fera jamais plus de 20 IA par veille

        try:
            ws_config = self.sh.worksheet(self.CONFIG_TAB_NAME)
            records = ws_config.get_all_records()

            configs = []
            seen_names = set()  # Pour détecter les doublons de noms

            for i, row in enumerate(records, start=2):
                is_active = str(row.get("Active", "")).upper() == "TRUE"

                if is_active:
                    raw_name = str(row.get("Nom Veille", "")).strip()

                    # --- SÉCURITÉ 1 : Protection du nom système ---
                    if raw_name == self.CONFIG_TAB_NAME:
                        logger.error(f"⛔ Ligne {i} : Nom '{self.CONFIG_TAB_NAME}' interdit.")
                        continue

                    # --- SÉCURITÉ 2 : Nettoyage du nom (Google compatible) ---
                    clean_name = sanitize_tab_name(raw_name)
                    if not clean_name:
                        logger.warning(f"⚠️ Ligne {i} : Nom vide ou invalide.")
                        continue

                    # --- SÉCURITÉ 3 : Détection de doublons ---
                    if clean_name in seen_names:
                        logger.warning(
                            f"⚠️ Ligne {i} : L'onglet '{clean_name}' est déjà défini plus haut. Ligne ignorée.")
                        continue
                    seen_names.add(clean_name)

                    # --- SÉCURITÉ 4 : Plafond de performance (Hard Cap) ---
                    user_max = int(row.get("Nb Max", 10))
                    # On prend le plus petit entre ce que veut l'user et notre limite technique
                    safe_max = min(user_max, max_allowed_results)

                    if user_max > max_allowed_results:
                        logger.warning(
                            f"ℹ️ Ligne {i} : Demande de {user_max} articles réduite à {max_allowed_results} pour performance.")

                    days = row.get("Jours Récents")
                    days_val = int(days) if days and str(days).isdigit() else None

                    try:
                        config = SearchConfig(
                            active=True,
                            name=clean_name,  # On utilise le nom nettoyé
                            query=row.get("Requête"),
                            max_results=safe_max,  # On utilise le max sécurisé
                            days_back=days_val,
                            start_date=str(row.get("Date Début", "")).strip() or None,
                            end_date=str(row.get("Date Fin", "")).strip() or None
                        )
                        configs.append(config)
                    except Exception as e:
                        logger.error(f"❌ Erreur format ligne {i} : {e}")
                        continue

            return configs

        except Exception as e:
            logger.error(f"Erreur lecture config : {e}")
            return []

    def set_target_tab(self, tab_name: str):
        """
        Change l'onglet de travail. Le crée s'il n'existe pas.
        """
        try:
            self.worksheet = self.sh.worksheet(tab_name)
        except gspread.WorksheetNotFound:
            logger.info(f"✨ Création de l'onglet '{tab_name}'...")
            self.worksheet = self.sh.add_worksheet(title=tab_name, rows=100, cols=10)
            # On initialise les en-têtes
            self.worksheet.append_row(["PMID", "Date", "Titre", "Score", "Résumé", "Points Clés", "URL"])
            self.update_layout()  # On applique le style immédiatement

    def get_existing_pmids(self) -> Set[str]:
        if not self.worksheet: return set()
        try:
            ids = self.worksheet.col_values(1)
            return set(ids[1:]) if ids else set()
        except Exception as e:
            logger.error(f"Erreur lecture PMIDs : {e}")
            return set()

    def save_article(self, enriched: EnrichedArticle):
        if not self.worksheet: return

        summary_text = enriched.summary.summary_text if enriched.summary else ""
        score = enriched.summary.relevance_score if enriched.summary else 0
        key_points = "\n".join(enriched.summary.key_points) if enriched.summary else ""

        row = [
            enriched.article.pubmed_id,
            str(enriched.article.publication_date),
            enriched.article.title,
            score,
            summary_text,
            key_points,
            enriched.article.url
        ]
        try:
            self.worksheet.append_row(row)
            logger.debug(f"💾 Sauvegardé : {enriched.article.pubmed_id}")
        except Exception as e:
            logger.error(f"Erreur écriture : {e}")

    def update_layout(self):
        """
        Applique le formatage via l'API BatchUpdate (Méthode Robuste).
        """
        if not self.worksheet: return

        try:
            print(f"   [Style] Application du formatage sur '{self.worksheet.title}'...")

            # On a besoin de l'ID unique de l'onglet pour l'API
            sheet_id = self.worksheet.id

            # Définition des largeurs de colonnes (Pixels)
            # Index 0 = Colonne A
            column_widths = [
                (0, 80),  # A: PMID
                (1, 90),  # B: Date
                (2, 250),  # C: Titre
                (3, 50),  # D: Score
                (4, 500),  # E: Résumé
                (5, 300),  # F: Points Clés
                (6, 200)  # G: URL
            ]

            requests = []

            # 1. FIGER LA LIGNE 1 (Freeze)
            requests.append({
                "updateSheetProperties": {
                    "properties": {
                        "sheetId": sheet_id,
                        "gridProperties": {"frozenRowCount": 1}
                    },
                    "fields": "gridProperties.frozenRowCount"
                }
            })

            # 2. LARGEUR DES COLONNES
            for col_index, width in column_widths:
                requests.append({
                    "updateDimensionProperties": {
                        "range": {
                            "sheetId": sheet_id,
                            "dimension": "COLUMNS",
                            "startIndex": col_index,
                            "endIndex": col_index + 1
                        },
                        "properties": {"pixelSize": width},
                        "fields": "pixelSize"
                    }
                })

            # 3. STYLE EN-TÊTE (Gras, Centré, Gris)
            requests.append({
                "repeatCell": {
                    "range": {
                        "sheetId": sheet_id,
                        "startRowIndex": 0,
                        "endRowIndex": 1
                    },
                    "cell": {
                        "userEnteredFormat": {
                            "backgroundColor": {"red": 0.9, "green": 0.9, "blue": 0.9},
                            "horizontalAlignment": "CENTER",
                            "textFormat": {"bold": True}
                        }
                    },
                    "fields": "userEnteredFormat(backgroundColor,textFormat,horizontalAlignment)"
                }
            })

            # 4. STYLE CORPS (Wrap Text + Alignement Haut)
            # On applique de la ligne 1 (la 2ème visuellement) jusqu'à la fin
            requests.append({
                "repeatCell": {
                    "range": {
                        "sheetId": sheet_id,
                        "startRowIndex": 1
                    },
                    "cell": {
                        "userEnteredFormat": {
                            "wrapStrategy": "WRAP",
                            "verticalAlignment": "TOP"
                        }
                    },
                    "fields": "userEnteredFormat(wrapStrategy,verticalAlignment)"
                }
            })

            # Envoi de la grosse commande à Google
            self.sh.batch_update({"requests": requests})
            print("   [Style] Mise en forme terminée avec succès.")

        except Exception as e:
            # On affiche l'erreur détaillée pour comprendre si l'API Google refuse
            print(f"⚠️ Attention: Le formatage a échoué : {e}")


def sanitize_tab_name(name: str) -> str:
    """
    Nettoie un nom pour qu'il soit accepté par Google Sheets.
    Règles : Pas de characters spéciaux interdits, max 31 chars.
    """
    # 1. Remplacer les caractères interdits par un underscore
    # Interdits : [ ] * ? : / \
    clean_name = re.sub(r"[\[\]*?:\\/]", "_", name)

    # 2. Couper à 31 caractères (limite stricte Google)
    clean_name = clean_name[:31]

    return clean_name