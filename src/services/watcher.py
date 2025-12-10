from src.infra.pubmed_client import PubMedClient
from src.infra.ai_client import OllamaClient
from src.infra.gsheets_repo import GoogleSheetsRepository
from src.domain.models import EnrichedArticle
from src.logger import logger


class SanitaryWatcher:
    def __init__(self):
        logger.info("--- Initialisation des services ---")
        self.pubmed_client = PubMedClient()
        self.ai_client = OllamaClient()
        self.repo = GoogleSheetsRepository()

    def run_all_watches(self):
        """
        Lit la configuration et exécute toutes les veilles actives.
        """
        # 1. Lecture de la config centralisée
        configs = self.repo.get_configs()

        if not configs:
            logger.warning("Aucune configuration active trouvée. Fin du programme.")
            return

        logger.info(f"🚀 Démarrage de la séquence : {len(configs)} veilles à traiter.")

        # 2. Boucle sur chaque ligne de configuration
        for config in configs:
            try:
                self._process_single_watch(config)
            except Exception as e:
                # ADR-009 : On log l'erreur mais on continue à la veille suivante
                logger.error(f"❌ CRASH CRITIQUE sur la veille '{config.name}' : {e}")
                continue

    def _process_single_watch(self, config):
        """
        Traite une seule veille de A à Z.
        """
        logger.info(f"\n=== Traitement Veille : {config.name} ===")

        # A. Préparer l'onglet Google Sheet
        self.repo.set_target_tab(config.name)

        # B. Récupérer le cache (IDs déjà traités dans CET onglet)
        existing_ids = self.repo.get_existing_pmids()

        # C. Recherche PubMed (Hybride : Query + Dates)
        found_ids = self.pubmed_client.search_article_ids(
            query=config.query,
            max_results=config.max_results,
            reldate=config.days_back,
            mindate=config.start_date,
            maxdate=config.end_date
        )

        # D. Filtrage
        new_ids = [pmid for pmid in found_ids if pmid not in existing_ids]

        if not new_ids:
            logger.info(f"✅ {config.name} : Rien de nouveau.")
            return

        logger.info(f"⚡ {len(new_ids)} nouveaux articles à traiter.")

        # E. Récupération & Analyse
        articles = self.pubmed_client.fetch_article_details(new_ids)

        count = 0
        for art in articles:
            try:
                summary = self.ai_client.summarize_article(art)
                enriched = EnrichedArticle(article=art, summary=summary)
                self.repo.save_article(enriched)
                count += 1
            except Exception as e:
                logger.error(f"Erreur article {art.pubmed_id} : {e}")

        # F. Formatage final
        if count > 0:
            self.repo.update_layout()
            logger.info(f"✅ Terminé pour {config.name} : {count} articles ajoutés.")