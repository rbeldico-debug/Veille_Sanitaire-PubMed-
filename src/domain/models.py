from datetime import date
from typing import Optional, List
from pydantic import BaseModel, Field

# Principe : Définir des objets clairs qui représentent votre métier.

class PubMedArticle(BaseModel):
    """
    Représente un article brut récupéré depuis PubMed.
    ADR-002 : Validation stricte des types.
    """
    pubmed_id: str = Field(..., description="L'identifiant unique PMID")
    title: str
    abstract: Optional[str] = None
    publication_date: Optional[date] = None
    authors: List[str] = Field(default_factory=list)
    url: str

class ArticleSummary(BaseModel):
    """
    Représente le résultat de l'analyse par l'IA.
    """
    summary_text: str
    relevance_score: int = Field(..., ge=0, le=10, description="Pertinence de 0 à 10")
    key_points: List[str] = Field(default_factory=list)

class EnrichedArticle(BaseModel):
    """
    Agrégation de l'article et de son analyse.
    C'est cet objet qui transitera vers Google Sheets.
    """
    article: PubMedArticle
    summary: Optional[ArticleSummary] = None
    processed_at: date = Field(default_factory=date.today)