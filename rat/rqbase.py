
import time
from rq import Queue
from redis import Redis

RC = Redis()
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


# #@rq_sync('t2', RQREDIS)
# def async(q, func, *args, **kwargs):
#     print('start async run')
#     job = q.enqueue(func, *args, **kwargs)
#     print(job.id)
#     while not job.is_finished:
#         time.sleep(1)
#     print('done rq run')
#     return job.result
