import gmpy2
import os

def read_number(filename):
    with open(filename, 'r') as f:
        return gmpy2.mpz(f.read().strip())

def main():
    base_path = r'C:\Users\Administrator\Documents\fft\\'
    num1_path = os.path.join(base_path, 'num1.txt')
    num2_path = os.path.join(base_path, 'num2.txt')
    result_path = os.path.join(base_path, 'result.txt')

    print(f"Reading {num1_path}...")
    num1 = read_number(num1_path)
    print(f"Reading {num2_path}...")
    num2 = read_number(num2_path)
    print(f"Reading {result_path}...")
    result = read_number(result_path)

    print("Calculating num1 * num2 with gmpy2...")
    expected = num1 * num2

    if expected == result:
        print("✅ Success: The result is correct!")
    else:
        print("❌ Error: The result does NOT match!")

if __name__ == '__main__':
    main()
