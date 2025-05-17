# propose only relevant files when imported
from . import analysis
from . import scoring
from . import preproc
from . import vizu

#Puts the emphasis on most relevant functions
from .preproc import cleanDic, NaCount
from .scoring import compute_scores

#For the autocompletion
__all__ = ["analysis", "scoring", "preproc", "vizu"]