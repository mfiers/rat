
import os
import time

from rq import Queue
from redis import Redis


_REDIS_ARGS = dict(
    host = os.environ.get('RQHOST'),
    port = int(os.environ.get('RQPORT', 6379)),
    password = os.environ.get('RQPASS'),
    )


RC = Redis(**_REDIS_ARGS)
QS = {}


def get_queue(name):
    if not name in QS:
        QS[name] = Queue(connection=RC, name=name)

    return QS[name]

def get_rq_redis_connection():
    return RC


def syncrun(func):
    """ decorator - to be used on top of rq's @job
        adds a `.sync` function to submit the job to
        rq, and then wait for completion 
    """
    def syncrunner(*args, **kwargs):
        import time
        job = func.delay(*args, **kwargs)
        while not job.is_finished and not job.is_failed:
            time.sleep(0.5)

        if job.is_failed:
            raise Exception('failed')

            
        print(dir(job))
        print('HEYYY', func)
        return job.result
        
    func.sync = syncrunner
    return func


