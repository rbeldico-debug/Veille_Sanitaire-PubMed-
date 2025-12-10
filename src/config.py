import os
from pathlib import Path
from dotenv import load_dotenv
from pydantic_settings import BaseSettings, SettingsConfigDict

# 1. Calcul du chemin et chargement explicite
PROJECT_DIR = Path(__file__).resolve().parent.parent
ENV_FILE_PATH = PROJECT_DIR / ".env"
load_dotenv(dotenv_path=ENV_FILE_PATH, override=True)

class Settings(BaseSettings):
    """
    Configuration de l'application.
    """
    PUBMED_EMAIL: str
    PUBMED_TOOL: str = "VeilleSanitaireApp"

    OLLAMA_BASE_URL: str = "http://localhost:11434"
    OLLAMA_MODEL: str = "gpt-oss:20b"

    GOOGLE_SHEET_NAME: str = "Veille_Sanitaire_Data"
    GOOGLE_CREDENTIALS_FILE: str = "google_credentials.json" # Le nom du fichier JSON
    
    # Configuration Pydantic V2
    model_config = SettingsConfigDict(
        env_file=str(ENV_FILE_PATH),
        env_file_encoding='utf-8',
        extra='ignore'
    )

settings = Settings()