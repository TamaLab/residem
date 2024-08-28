import concurrent.futures
import os
from functools import wraps

import time

def make_parallel(func):
    """
    Decorator used to decorate any function which needs to be parallelized.
    After the input of the function should be a list in which each element is an instance of input for the normal function.
    You can also pass in keyword arguments separately.
    :param func: function
        The instance of the function that needs to be parallelized.
    :return: function
    """
    @wraps(func)
    def wrapper(lst):
        """
        :param lst:
            The inputs of the function in a list.
        :return:
        """
        result = []
        while True:
            try:
                number_of_threads_multiple = 2  # You can change this multiple according to your requirement
                number_of_workers = int(os.cpu_count() * number_of_threads_multiple)
                if len(lst) < number_of_workers:
                    number_of_workers = len(lst)

                if number_of_workers > 1:
                    with concurrent.futures.ThreadPoolExecutor(max_workers=number_of_workers) as executor:
                        bag = {executor.submit(func, i): i for i in lst}
                        for future in concurrent.futures.as_completed(bag):
                            result.append(future.result())
                elif number_of_workers == 1:
                    result = [func(lst[0])]
                else:
                    result = []

                break  # Break the loop if the code executes without a RuntimeError
            except RuntimeError as e:
                print(f"Error occurred: {e}")
                print("Retrying after 10 seconds...")
                time.sleep(10)  # Wait for 5 seconds before retrying

        return result

    return wrapper
