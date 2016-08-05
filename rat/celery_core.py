


import os
from celery import Celery


_APP = None

def get_celery_app():
    if _APP is None:
        global _APP
        broker = os.environ.get('CELERY_BROKER', 'redis://:muffins1@r10n1:6379/0')
        backend = os.environ.get('CELERY_BACKEND', 'redis://:muffins1@r10n1:6379/0')
        concurrency = int(os.environ.get('CELERYD_CONCURRENCY', 11))

        broker_transport_options = {'fanout_patterns': True,
                                    'fanout_prefix': True}

        _APP = Celery('tasks', broker=broker, backend=backend)
        _APP.conf.CELERYD_POOL_RESTARTS = True
        _APP.conf.CELERYD_CONCURRENCY = concurrency
        _APP.conf.CELERY_ACKS_LATE = True
    return _APP
