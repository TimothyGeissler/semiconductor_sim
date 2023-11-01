def memo(n):
    saved = [0] * (n + 1)
    return fib(n, saved)

def fib(n, saved):
    if (n == 0 or n == 1):
        return 1
    if (saved[n] == 0):
        saved[n] = fib(n - 1, saved) + fib(n - 2, saved)
    return saved[n]

def fib_slo(n):
    if (n == 0 or n == 1):
        return 1
    return fib_slo(n - 1) + fib_slo(n - 2)

# bottom up DP
def dpfib(n):
    a = 1
    b = 1
    for i in range(1, n):
       c = a + b
       a = b
       b = c
    return c 
    
print(memo(100))
print(dpfib(100))