

from functools import partial
import logging
import multiprocessing as mp
import re
import sys
import socket
import time

from rat.rqbase import get_rq_redis_connection

from rq import Connection, Worker
from rq import use_connection

lg = logging.getLogger(__name__)


class NeedMoreSpaceWorker(Worker):

    def execute_job(self, job, queue):

        self.set_state('busy')
        
        JID = job.id
        qname = queue.name
        qthread = int(qname[1:])

        #add qthread tokes to the local queue
        self.lqueue.extend([JID] * qthread)

        def maxindex(lst, val):
            for i, l in enumerate(lst[::-1]):
                if l == val: return len(lst) - i - 1

        # now wait until the hightest index token < max no slots:
        while maxindex(self.lqueue, JID) > self.maxthreads:
		            time.sleep(0.5)

        print('GOGO', qname, self.maxthreads, '# midx',
              maxindex(self.lqueue, JID), 'qlen', len(self.lqueue))
        print('   -', JID)

        self.set_state('busy')

        self.fork_work_horse(job, queue)
        self.monitor_work_horse(job)

        #remove all tokens from queue
        try:
            while True:
                self.lqueue.remove(JID)
        except ValueError:
            pass
            
        print("DONE", qname, JID)
        self.set_state('idle')
										

def SUPERWORKER(i, maxthreads, rqueues, lqueue):
    hostname = socket.gethostname()
    workername = '{:8s}.{:2d}.{:2d}'.format(hostname, maxthreads, i).replace(' ', '.')
    with Connection(get_rq_redis_connection() ):
        w = NeedMoreSpaceWorker(rqueues, name=workername)
        w.lqueue = lqueue
        w.maxthreads = maxthreads
        w.work()

def dispatch():
    maxthreads=int(sys.argv[1])
    QUEUES = ['t{}'.format(i) for i in range(1, maxthreads+1)]

    pool = mp.Pool(maxthreads)
    manager = mp.Manager()
    lqueue = manager.list() #shared between states
    swp = partial(SUPERWORKER, maxthreads=maxthreads, rqueues=QUEUES, lqueue=lqueue)
    pool.map(swp, range(1, maxthreads+1))
    pool.close()
    pool.join()

                                                    
