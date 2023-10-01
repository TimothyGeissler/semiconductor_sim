def is_palindrome(s):
    return s == s[::-1]

def min_palindromes_to_cover(s):
    n = len(s)
    dp = [float('inf')] * n  # Initialize an array to store minimum counts
    dp[0] = 1  # The first character is always a palindrome

    for i in range(1, n):
        for j in range(i, -1, -1):
            if is_palindrome(s[j:i+1]):
                if j == 0:
                    dp[i] = 1 # 0 -> i can be covered in 1 palindrome
                else:
                    dp[i] = min(dp[i], dp[j - 1] + 1) #min between current value and 1 + palin count at the beginning of current palindrome

    print(dp)
    return dp[n - 1]

# Example usage:
input_string = "BUBBASEESABANANA"
result = min_palindromes_to_cover(input_string2 )
print(result)