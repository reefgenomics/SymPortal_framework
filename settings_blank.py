import os

dbPath = os.path.join(os.path.dirname(__file__), 'db.sqlite3')

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3', # Add 'postgresql_psycopg2', 'mysql', 'sqlite3' or 'oracle'.
        'NAME': '{}'.format(dbPath),                      # Or path to database file if using sqlite3.
        'OPTIONS': {'timeout':20}
    }
}



# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = False

CACHES = {
    'default': {
        'BACKEND': 'django.core.cache.backends.dummy.DummyCache',
    }
}


INSTALLED_APPS = (
    'dbApp',
    )

SECRET_KEY = ''
