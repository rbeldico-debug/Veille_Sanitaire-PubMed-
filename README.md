# Veille Sanitaire Automatisée 🦠

Ce projet est un outil de veille scientifique automatisé.
Il surveille **PubMed** pour des mots-clés spécifiques, analyse les nouveaux articles via une **IA locale (Ollama)** pour en extraire l'essentiel, et publie les résultats formatés dans **Google Sheets**.

## 🚀 Fonctionnalités

*   **100% Configurable** : Les veilles se définissent directement dans un Google Sheet (pas besoin de toucher au code).
*   **Intelligent** : Résume les articles et attribue un score de pertinence via IA.
*   **Incrémental** : Ne traite jamais deux fois le même article (détection de doublons).
*   **Hybride** : Recherche par mots-clés sémantiques ET par dates.
*   **Robuste** : Logs rotatifs et gestion d'erreurs (une veille plantée ne bloque pas les autres).

## 🛠 Prérequis Techniques

1.  **Python 3.10+**
2.  **Ollama** installé localement (`ollama serve`) avec le modèle configuré (ex: `gpt-oss:20b`).
3.  Un compte de service **Google Cloud** (fichier `google_credentials.json` à la racine).

## ⚙️ Configuration (Utilisateur)

Tout se passe dans le Google Sheet maître, dans l'onglet **`_ADMIN_CONFIG`**.
Chaque ligne correspond à une veille. Le script créera automatiquement les onglets de résultats.

| Colonne | Description | Exemple |
| :--- | :--- | :--- |
| **Active** | `TRUE` ou `FALSE`. Active/Désactive la veille. | `TRUE` |
| **Nom Veille** | Nom de l'onglet qui sera créé/rempli. | `Grippe Aviaire` |
| **Requête** | Mots-clés PubMed (Syntaxe standard). | `H5N1 AND France` |
| **Nb Max** | Nombre max d'articles à analyser par exécution. | `10` |
| **Jours Récents** | Nombre de jours en arrière (Prioritaire). | `30` |
| **Date Début** | Format YYYY/MM/DD (Si Jours Récents est vide). | `2023/01/01` |
| **Date Fin** | Format YYYY/MM/DD. | `2023/12/31` |

## 💻 Installation (Développeur)

1.  Cloner le dépôt :
    ```bash
    git clone https://github.com/VOTRE_USER/Veille_Sanitaire.git
    ```
2.  Créer l'environnement virtuel et installer les dépendances :
    ```bash
    python -m venv .venv
    .venv\Scripts\activate
    pip install -r requirements.txt
    ```
3.  Configurer les secrets :
    *   Créer un fichier `.env` à la racine (voir `src/config.py`).
    *   Placer `google_credentials.json` à la racine.

## ▶️ Exécution

*   **Manuelle** : Double-cliquer sur `run.bat`.
*   **Automatique** : Configurer le *Planificateur de tâches Windows* pour lancer `run.bat` quotidiennement.

## 📁 Logs

L'exécution est tracée dans le fichier `execution.log` à la racine.