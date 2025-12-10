from typing import List
from src.domain.models import EnrichedArticle
from src.infra.pubmed_client import PubMedClient
from src.infra.ai_client import OllamaClient
from src.infra.gsheets_repo import GoogleSheetsRepository


class SanitaryWatcher:
    """
    Orchestrateur principal de la veille sanitaire.
    Coordonne PubMed, l'IA et Google Sheets.
    """

    def __init__(self):
        # Initialisation des 3 dépendances d'infrastructure
        print("Initialisation des services...")
        self.pubmed_client = PubMedClient()
        self.ai_client = OllamaClient()
        self.repo = GoogleSheetsRepository()

    def run_process(self, query: str, max_results: int = 10):
        """
        Exécute le workflow complet de veille.
        """
        print(f"\n=== Lancement de la veille pour '{query}' ===")

        # ÉTAPE 1 : Récupérer ce qu'on connait déjà (Cache)
        print("1. Lecture de la base de données existante...")
        existing_ids = self.repo.get_existing_pmids()
        print(f"   -> {len(existing_ids)} articles déjà stockés.")

        # ÉTAPE 2 : Chercher les IDs récents sur PubMed
        print(f"2. Recherche des derniers articles sur PubMed (max {max_results})...")
        found_ids = self.pubmed_client.search_article_ids(query, max_results)

        # ÉTAPE 3 : Filtrage (Logique de dédoublonnement)
        # On garde l'ID seulement s'il N'EST PAS dans existing_ids
        new_ids = [pmid for pmid in found_ids if pmid not in existing_ids]

        if not new_ids:
            print("\n✅ Aucune nouveauté détectée. Tout est à jour.")
            return

        print(f"   -> {len(new_ids)} nouveaux articles identifiés à traiter.")

        # ÉTAPE 4 : Récupération des détails (XML) pour les nouveaux seulement
        articles_to_process = self.pubmed_client.fetch_article_details(new_ids)

        # ÉTAPE 5 : Boucle de traitement (IA + Sauvegarde)
        print(f"\n3. Analyse IA et Sauvegarde ({len(articles_to_process)} articles)...")

        count_success = 0
        for i, article in enumerate(articles_to_process, 1):
            print(f"\n--- Traitement {i}/{len(articles_to_process)} : PMID {article.pubmed_id} ---")
            try:
                # A. Résumé IA
                summary = self.ai_client.summarize_article(article)

                # B. Création de l'objet enrichi
                enriched = EnrichedArticle(article=article, summary=summary)

                # C. Sauvegarde immédiate (comme ça si ça plante au 10ème, les 9 premiers sont sauvés)
                self.repo.save_article(enriched)
                count_success += 1

            except Exception as e:
                print(f"❌ Erreur sur l'article {article.pubmed_id} : {e}")

        if count_success > 0:
            self.repo.update_layout()

        print(f"\n=== Terminé : {count_success} articles ajoutés au Google Sheet. ===")