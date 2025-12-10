import sys
from src.services.watcher import SanitaryWatcher


def main():
    # Configuration de la recherche
    # Vous pourrez plus tard passer ces paramètres via la ligne de commande
    QUERY = "Psychiatry"
    MAX_RESULTS_TO_CHECK = 5  # On vérifie les 5 derniers articles sortis

    try:
        # Instanciation du service
        watcher = SanitaryWatcher()

        # Lancement du processus
        watcher.run_process(query=QUERY, max_results=MAX_RESULTS_TO_CHECK)

    except KeyboardInterrupt:
        print("\nArrêt manuel par l'utilisateur.")
    except Exception as e:
        print(f"\nCRASH GLOBAL : {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()