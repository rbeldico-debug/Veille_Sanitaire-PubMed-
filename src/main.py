import sys
from src.services.watcher import SanitaryWatcher
from src.logger import logger

def main():
    try:
        watcher = SanitaryWatcher()
        watcher.run_all_watches()
    except Exception as e:
        logger.critical(f"Erreur fatale au lancement : {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()