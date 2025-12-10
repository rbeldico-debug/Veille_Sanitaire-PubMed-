import gspread
from typing import Set
from src.config import settings, PROJECT_DIR
from src.domain.models import EnrichedArticle


class GoogleSheetsRepository:
    """
    Gère la persistance des données dans Google Sheets.
    """

    def __init__(self):
        # Chemin absolu vers le fichier JSON de credentials
        creds_path = PROJECT_DIR / settings.GOOGLE_CREDENTIALS_FILE

        if not creds_path.exists():
            raise FileNotFoundError(f"Le fichier de crédentials Google est introuvable : {creds_path}")

        # Connexion à Google
        self.gc = gspread.service_account(filename=str(creds_path))

        # Ouverture du fichier (Spreadsheet)
        try:
            self.sh = self.gc.open(settings.GOOGLE_SHEET_NAME)
            # On suppose qu'on écrit dans la première feuille (Sheet1)
            self.worksheet = self.sh.sheet1
        except gspread.SpreadsheetNotFound:
            raise ValueError(
                f"Impossible de trouver le Google Sheet nommé : '{settings.GOOGLE_SHEET_NAME}'. Vérifiez le nom et le partage avec le robot.")

    def get_existing_pmids(self) -> Set[str]:
        """
        Récupère tous les IDs (colonne A) pour le dédoublonnement.
        """
        try:
            # Récupère toute la colonne 1 (PMID)
            ids = self.worksheet.col_values(1)
            if ids:
                # On retire le header de la ligne 1 pour ne garder que les IDs
                return set(ids[1:])
            return set()
        except Exception as e:
            print(f"Erreur lecture GSheets : {e}")
            return set()

    def save_article(self, enriched: EnrichedArticle):
        """
        Ajoute une ligne à la fin du tableau.
        """
        summary_text = enriched.summary.summary_text if enriched.summary else ""
        score = enriched.summary.relevance_score if enriched.summary else 0
        key_points = "\n".join(enriched.summary.key_points) if enriched.summary else ""

        # Ordre des colonnes : PMID | Date | Titre | Score | Résumé | Points Clés | URL
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
            print(f"   [GSheets] Sauvegardé : {enriched.article.pubmed_id}")
        except Exception as e:
            print(f"Erreur écriture GSheets : {e}")

    def update_layout(self):
        """
        Applique un formatage visuel pour rendre le tableau lisible.
        """
        try:
            print("   [Style] Mise en forme du tableau en cours...")

            # 1. Figer la ligne d'en-tête
            self.worksheet.freeze(rows=1)

            # 2. Largeur des colonnes (Nom exact : set_column_width)
            # Attention : gspread utilise des indices commençant à 1 pour A
            self.worksheet.set_column_width(1, 100)  # A: PMID
            self.worksheet.set_column_width(2, 100)  # B: Date
            self.worksheet.set_column_width(3, 250)  # C: Titre
            self.worksheet.set_column_width(4, 50)  # D: Score
            self.worksheet.set_column_width(5, 500)  # E: Résumé (Large)
            self.worksheet.set_column_width(6, 300)  # F: Points Clés
            self.worksheet.set_column_width(7, 150)  # G: URL

            # 3. Application des styles (Gras pour titre, Wrap pour le texte)
            # On définit le style pour tout le tableau d'un coup

            # En-tête : Gras + Centré + Fond gris clair
            self.worksheet.format("A1:G1", {
                "textFormat": {"bold": True},
                "horizontalAlignment": "CENTER",
                "backgroundColor": {"red": 0.9, "green": 0.9, "blue": 0.9}
            })

            # Corps du texte (Colonnes C à F) : Retour à la ligne (WRAP) + Alignement Haut
            # On applique sur une plage assez grande (ex: ligne 2 à 1000)
            self.worksheet.format(f"C2:F{self.worksheet.row_count}", {
                "wrapStrategy": "WRAP",
                "verticalAlignment": "TOP"
            })

            print("   [Style] Mise en forme terminée.")

        except Exception as e:
            # On affiche l'erreur mais on ne bloque pas le programme
            print(f"Attention: Le formatage visuel a échoué : {e}")