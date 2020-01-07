def in_notebook():
    if hasattr(__builtins__,'__IPYTHON__'):
        return True #globals()[notebook_mod] = __import__(notebook_mod)
    else:
        return False #globals()[console_mod] = __import__(console_mod)

if in_notebook():
    from tqdm import tqdm_notebook as tqdm
else:
    from tqdm import tqdm as tqdm
