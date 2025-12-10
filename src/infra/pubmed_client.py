from typing import List, Optional
from datetime import date
from Bio import Entrez
from src.config import settings
from src.domain.models import PubMedArticle
from src.logger import logger


class PubMedClient:
    """
    Gère les interactions avec l'API PubMed via BioPython.
    Sépare la recherche d'IDs (esearch) de la récupération de détails (efetch)
    pour permettre le filtrage des doublons.
    """

    def __init__(self):
        # Configuration obligatoire pour l'API NCBI (ADR-003)
        Entrez.email = settings.PUBMED_EMAIL
        Entrez.tool = settings.PUBMED_TOOL

    def search_article_ids(self, query: str, max_results: int = 5,
                           reldate: Optional[int] = None,
                           mindate: Optional[str] = None,
                           maxdate: Optional[str] = None) -> List[str]:
        """
        Recherche hybride : Mots-clés + Dates (ADR-008).
        Priorité : reldate > mindate/maxdate.
        """
        try:
            params = {
                "db": "pubmed",
                "term": query,
                "retmax": max_results,
                "sort": "date",
                "retmode": "xml"
            }

            # Gestion des dates (ADR-008)
            date_info = ""
            if reldate:
                params["reldate"] = reldate
                params["datetype"] = "pdat"  # Date de publication
                date_info = f"(Derniers {reldate} jours)"
            elif mindate and maxdate:
                params["mindate"] = mindate
                params["maxdate"] = maxdate
                params["datetype"] = "pdat"
                date_info = f"({mindate} -> {maxdate})"

            logger.info(f"🔎 Recherche PubMed : '{query}' {date_info}")

            # Appel API avec unpack des arguments (**params)
            search_handle = Entrez.esearch(**params)
            search_results = Entrez.read(search_handle)
            search_handle.close()

            id_list = search_results["IdList"]
            logger.info(f"   -> {len(id_list)} IDs trouvés.")
            return id_list

        except Exception as e:
            logger.error(f"❌ Erreur lors de la recherche IDs : {e}")
            return []

    def fetch_article_details(self, id_list: List[str]) -> List[PubMedArticle]:
        """
        Équivalent des Noeuds n8n 'HTTP Request (efetch)' + 'XML Parse'.
        Prend une liste d'IDs (déjà filtrée par l'orchestrateur) et récupère les détails.
        """
        if not id_list:
            return []

        try:
            print(f"Téléchargement des détails pour {len(id_list)} articles...")

            # Conversion de la liste Python en chaîne séparée par des virgules pour l'API
            ids_str = ",".join(id_list)

            fetch_handle = Entrez.efetch(
                db="pubmed",
                id=ids_str,
                retmode="xml"
            )
            # Entrez.read parse le XML automatiquement (équivalent du noeud XML n8n)
            articles_data = Entrez.read(fetch_handle)
            fetch_handle.close()

            pubmed_articles = []
            if 'PubmedArticle' in articles_data:
                for article_xml in articles_data['PubmedArticle']:
                    domain_article = self._map_xml_to_domain(article_xml)
                    if domain_article:
                        pubmed_articles.append(domain_article)

            return pubmed_articles

        except Exception as e:
            print(f"Erreur lors du fetch details : {e}")
            return []

    def _map_xml_to_domain(self, xml_data) -> Optional[PubMedArticle]:
        """
        Convertit le dictionnaire brut de BioPython en PubMedArticle (Pydantic).
        """
        try:
            medline = xml_data['MedlineCitation']
            article_info = medline['Article']

            pmid = str(medline['PMID'])
            title = article_info.get('ArticleTitle', 'Titre inconnu')

            # Gestion de l'Abstract
            abstract_text = "Pas de résumé disponible."
            if 'Abstract' in article_info and 'AbstractText' in article_info['Abstract']:
                abstract_source = article_info['Abstract']['AbstractText']
                if isinstance(abstract_source, list):
                    # Parfois l'abstract est découpé en sections (Background, Methods...)
                    # On concatène tout proprement
                    abstract_text = " ".join(str(item) for item in abstract_source)
                else:
                    abstract_text = str(abstract_source)

            # Gestion des Auteurs
            authors = []
            if 'AuthorList' in article_info:
                for author in article_info['AuthorList']:
                    # Parfois les auteurs sont des collectifs sans ForeName/LastName
                    if 'LastName' in author and 'ForeName' in author:
                        authors.append(f"{author['ForeName']} {author['LastName']}")
                    elif 'LastName' in author:
                        authors.append(author['LastName'])
                    elif 'CollectiveName' in author:
                        authors.append(author['CollectiveName'])

            pub_date = self._extract_date(article_info)
            url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

            return PubMedArticle(
                pubmed_id=pmid,
                title=title,
                abstract=abstract_text,
                publication_date=pub_date,
                authors=authors,
                url=url
            )
        except Exception as e:
            # On log l'erreur mais on ne plante pas tout le processus pour un article malformé
            print(f"Erreur parsing article {xml_data.get('MedlineCitation', {}).get('PMID', '?')}: {e}")
            return None

    def _extract_date(self, article_info) -> Optional[date]:
        """Helper pour extraire une date."""
        try:
            if 'ArticleDate' in article_info and len(article_info['ArticleDate']) > 0:
                d = article_info['ArticleDate'][0]
                return date(int(d['Year']), int(d['Month']), int(d['Day']))

            journal_date = article_info.get('Journal', {}).get('JournalIssue', {}).get('PubDate', {})
            if 'Year' in journal_date:
                return date(int(journal_date['Year']), 1, 1)

            return None
        except Exception:
            return None