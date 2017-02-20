# -*- coding: utf-8 -*-
"""
Created on Mon Feb 06 16:04:07 2017

@author: afalaize
"""

from __future__ import division, absolute_import, print_function

import pp
import numpy as np
from multiprocessing import Pool, Queue
import time
from contextlib import closing


import threading



ncpus = 2
def mapfunc(func, argslist, deps=tuple(), modules=tuple()):
    """
    usage: 
        res = mapfunc(func, argslist, deps=tuple(), modules=tuple())
    """
    job_server = pp.Server(ncpus)
    print('{} cpus detected'.format(job_server.get_ncpus()))    
    jobs = list()
    for arg in argslist:
        jobs.append(job_server.submit(func, arg, deps, modules))
    results = list()
    for job in jobs:
        results.append(job())
    print('\nStats for parallelization:')
    job_server.print_stats()
    return results





def threadingmap(func, listOfArgs):
    
    q = Queue()
    
    def putTaskInQueue(args, q):
        q.put(func(args))
    
    tasks = list()
    
    for args in listOfArgs:
        tasks.append(threading.Thread(target=putTaskInQueue, 
                                      args=(args, q)))

    for task in tasks:
        task.start()
        
    for task in tasks:
        task.join()
    def output_generator():
        while not q.empty():
            yield q.get()
    return list(output_generator())
    
import concurrent


def concurentmap(func, listOfArgs, max_workers=16):
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        return list(executor.map(func, listOfArgs))

    
if __name__ == "__main__":        
    def add(t):
        return sum(t)
        
    args = [(i, 0) for i in range(100)]            
#    time_pp = time.time()
#    results = mapfunc(add, args)
#    print('pp took {}s'.format(time.time()-time_pp))
#   
 

    time_pool = time.time()
    res = threadingmap(add, args)
    print('threadingmap took {}s'.format(time.time()-time_pool))
    
    time_concurentmap = time.time()
    res = concurentmap(add, args)
    print('concurentmap took {}s'.format(time.time()-time_concurentmap))
        