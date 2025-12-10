# Veille Sanitaire Automatisée

Ce projet surveille PubMed pour des mots clés spécifiques (ex: H5N1), analyse les nouveaux articles via une IA locale (Ollama), et ajoute un résumé structuré dans un Google Sheet.

## Prérequis

1.  **Python 3.10+**
2.  **Ollama** installé et tournant (`ollama serve`) avec le modèle mistral (`ollama pull mistral`).
3.  Un compte **Google Cloud** avec un Service Account (fichier JSON).

## Installation

1.  Cloner le repo.
2.  Créer l'environnement virtuel : `python -m venv .venv`
3.  Activer : `.venv\Scripts\activate`
4.  Installer les dépendances : `pip install -r requirements.txt`
5.  Placer le fichier `google_credentials.json` à la racine.
6.  Créer le fichier `.env` (voir `src/config.py` pour les variables requises).

## Utilisation

- **Lancement manuel** : Double-cliquer sur `run.bat`.
- **Automatisation** : Le script est prévu pour être lancé par le Planificateur de tâches Windows.

## Structure

- `src/domain` : Modèles de données.
- `src/infra` : Clients externes (PubMed, Google, Ollama).
- `src/services` : Logique de l'orchestrateur.