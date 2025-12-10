import logging
from logging.handlers import RotatingFileHandler
from src.config import PROJECT_DIR

# Chemin du fichier de log (à la racine)
LOG_FILE = PROJECT_DIR / "execution.log"


def setup_logger():
    """
    Configure un logger qui écrit dans un fichier ET dans la console.
    Rotation : Max 1Mo par fichier, garde 3 archives.
    """
    logger = logging.getLogger("VeilleSanitaire")
    logger.setLevel(logging.INFO)

    # Éviter les doublons si la fonction est appelée plusieurs fois
    if logger.handlers:
        return logger

    # 1. Handler Fichier (Rotatif)
    file_handler = RotatingFileHandler(
        LOG_FILE, maxBytes=1_000_000, backupCount=3, encoding="utf-8"
    )
    file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(file_formatter)

    # 2. Handler Console (Pour voir ce qui se passe en direct)
    console_handler = logging.StreamHandler()
    console_formatter = logging.Formatter('%(message)s')  # Plus simple pour la console
    console_handler.setFormatter(console_formatter)

    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger


# Instance globale
logger = setup_logger()