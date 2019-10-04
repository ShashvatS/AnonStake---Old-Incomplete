import secrets
from mpmath import mp, mpf, log, exp

MIMC_ROUNDS = 162
prime_r = 52435875175126190479447740508185965837690552500527637822603658699938581184513

cipherlist = ["hash", "addr", "sn", "tsn", "seed", "priority"] + ["binom" + str(i) for i in range(60)]

file = open("constants2.txt", "w")

FORMAT_FIELDR_NAMED = "const FieldR {}(\"{}\");"
FORMAT_FIELDR = "FieldR(\"{}\")"

file.write(FORMAT_FIELDR_NAMED.format("hashconstant", str(secrets.randbits(256) % prime_r)) + "\n")

for cipher in cipherlist:
    constants = [0] + [secrets.randbits(256) % prime_r for i in range(MIMC_ROUNDS - 1)]
    data = "const FieldR {}[{}] = {{ ".format("PRF" + cipher, MIMC_ROUNDS)

    for i in range(MIMC_ROUNDS - 1):
        data += FORMAT_FIELDR.format(constants[i]) + ", "
    data += FORMAT_FIELDR.format(constants[-1])
    data += "};\n"

    file.write(data)

# a bit extra
mp.prec = 80

numcoins = mpf(2**60 - 1)
prob = mpf(2000) / numcoins
negl = mpf(1.0) * pow(mpf(0.5), 80)

print(numcoins, prob, negl)

def log_nck(n, k):
    if k == 0 or k == n:
        return log(mpf(1))

    nint = n
    kint = k
    n = mpf(n)
    k = mpf(k)
    ans = mpf(0)
    for i in range(nint, nint - kint, -1):
        ans += log(mpf(i))
    for i in range(kint, 0, -1):
        ans -= log(mpf(i))
    return ans

def calc_prob(choose_i, from_n):
    # probability = (nCk) * p**k * (1-p)**(n - k)
    log_ans = log_nck(from_n, choose_i) + mpf(choose_i) * log(prob) 
    log_ans += (mpf(from_n - choose_i)) * log(mpf(1) - prob)
    return exp(log_ans)

for i in range(0, 60):
    n = pow(2, i)
    cum_prob = calc_prob(0, n)
    k = 0

    values = []

    while cum_prob < negl:
        k += 1
        cum_prob += calc_prob(k, n)
    k1 = k

    while cum_prob < mpf(1.0) - negl and k < n:
        values.append((k, cum_prob))
        k += 1
        tmp = calc_prob(k, n)
        cum_prob += tmp
        # print(cum_prob, end=" ")
        
    k2 = k
    print(k1, k2 - 1)

    values = [(x[0], int(x[1] * (mpf(2) ** 80))) for x in values]

    data = ["probvalue" + str(i)] + [str(values[0][0])] + [str(i[1]) for i in values] + ["\n"]
    file.write("qqqArr[{}] = {};\n".format(i, len(values)))
    file.write("qqqqArr[{}] = {};\n".format(i, values[0][0]))
    file.write("qqArr[{}] = new FieldR[{}];\n".format(i, len(values)))
    for j in range(len(values)):
        a = "qqArr[{}][{}] = "
        a = a.format(str(i), str(j))
        a += FORMAT_FIELDR.format(data[j + 2]) + ';\n'
        file.write(a)

file.close()






