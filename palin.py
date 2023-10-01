def is_palindrome(s):
    return s == s[::-1]

def find_palindromic_substrings(s):
    palindromes = []
    n = len(s)

    for i in range(n):
        # Odd-length palindromes centered at s[i]
        l, r = i, i
        while l >= 0 and r < n and s[l] == s[r]:
            if is_palindrome(s[l:r+1]):
                palindromes.append(s[l:r+1])
            l -= 1
            r += 1

        # Even-length palindromes centered between s[i] and s[i+1]
        l, r = i, i + 1
        while l >= 0 and r < n and s[l] == s[r]:
            if is_palindrome(s[l:r+1]):
                palindromes.append(s[l:r+1])
            l -= 1
            r += 1

    return palindromes

# Example usage:
input_string = "BUBBASEESABANANA"
palindromes = find_palindromic_substrings(input_string)
print(palindromes)
