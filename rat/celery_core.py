
import os
from celery import Celery

_APP = None

def get_celery_app():
    global _APP
    if _APP is None:

        broker = os.environ.get('CELERY_BROKER', 'redis://localhost:6379/0')
        backend = os.environ.get('CELERY_BACKEND', 'redis://localhost:6379/0')

        #broker = os.environ.get('CELERY_BROKER', 'redis://:muffins1@localhost:6379/0')
        #backend = os.environ.get('CELERY_BACKEND', 'redis://:muffins1@localhost:6379/0')

        #broker = os.environ.get('CELERY_BROKER', 'redis://:muffins1@r10n1:6379/0')
        #backend = os.environ.get('CELERY_BACKEND', 'redis://:muffins1@r10n1:6379/0')
        concurrency = int(os.environ.get('CELERYD_CONCURRENCY', 11))

        broker_transport_options = {'fanout_patterns': True,
                                    'fanout_prefix': True}

        _APP = Celery('tasks', broker=broker, backend=backend)
        _APP.conf.worker_pool_restarts = True
        _APP.conf.worker_concurrency = concurrency
        _APP.conf.task_acks_late = True
    return _APP
