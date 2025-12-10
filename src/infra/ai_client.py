import json
from langchain_ollama import ChatOllama  # <-- NOUVEL IMPORT (V2)
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.output_parsers import StrOutputParser

from src.logger import logger
from src.config import settings
from src.domain.models import PubMedArticle, ArticleSummary


class OllamaClient:
    """
    Client pour interagir avec un LLM local via Ollama.
    """

    def __init__(self):
        # Utilisation de la nouvelle classe ChatOllama
        self.llm = ChatOllama(
            base_url=settings.OLLAMA_BASE_URL,
            model=settings.OLLAMA_MODEL,
            format="json",
            temperature=0
        )

    def summarize_article(self, article: PubMedArticle) -> ArticleSummary:
        print(f"   Analysé par IA ({settings.OLLAMA_MODEL}) : {article.title[:50]}...")
        logger.debug(f"🤖 Analyse IA en cours : {article.title[:30]}...")
        # CORRECTION CRITIQUE ICI :
        # Les accolades du JSON exemple doivent être doublées {{ }}
        # sinon LangChain pense que ce sont des variables à remplacer.
        system_prompt = """
        Tu es un expert en veille sanitaire et épidémiologique.
        Ton objectif est d'analyser des articles scientifiques pour des professionnels de santé.

        Analyse l'article fourni et retourne UNIQUEMENT un objet JSON valide respectant exactement cette structure :
        {{
            "summary_text": "Résumé concis en français de l'article (max 3 phrases).",
            "relevance_score": 8, 
            "key_points": ["Point clé 1", "Point clé 2", "Point clé 3"]
        }}

        Le "relevance_score" doit être un entier de 0 à 10 (10 = crucial pour la santé publique, 0 = hors sujet).
        Traduis tout en français.
        """

        # On passe les variables proprement via le template
        prompt = ChatPromptTemplate.from_messages([
            ("system", system_prompt),
            ("user", "Titre : {title}\nRésumé brut : {abstract}")
        ])

        chain = prompt | self.llm | StrOutputParser()

        try:
            # On injecte les variables ici
            json_response = chain.invoke({
                "title": article.title,
                "abstract": article.abstract or "Pas de résumé"
            })

            data = json.loads(json_response)

            return ArticleSummary(
                summary_text=data.get("summary_text", "Erreur de résumé"),
                relevance_score=data.get("relevance_score", 0),
                key_points=data.get("key_points", [])
            )

        except Exception as e:
            print(f"Erreur lors de l'analyse IA : {e}")
            return ArticleSummary(
                summary_text="Erreur technique lors de l'analyse IA.",
                relevance_score=0,
                key_points=[]
            )