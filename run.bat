@echo off
:: ==========================================
:: Script de lancement automatique Veille Sanitaire
:: ==========================================

:: 1. Se déplacer dans le dossier du projet (Chemin Absolu)
cd /d "C:\Users\G-i7\PycharmProjects\Veille_Sanitaire"

:: 2. Activer l'environnement virtuel (Le dossier .venv)
call .venv\Scripts\activate.bat

:: 3. Lancer le programme Python
:: On utilise -m pour lancer le module proprement
echo Lancement de la veille...
python -m src.main

:: 4. Attendre 30 secondes avant de fermer la fenêtre
:: (Permet de lire les erreurs s'il y en a au démarrage)
timeout /t 30