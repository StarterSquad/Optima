import sys
import os

try:
    import server
except:
    sys.path.insert(0, os.path.abspath(".."))

from server import _autoreload
from server import _twisted_wsgi

_autoreload.main(_twisted_wsgi.run)
