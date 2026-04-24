MODULO = 17
POWER = 1729

INTEGER  = 1729

step = 0
remainder_old = INTEGER
remainder_new = -1


while remainder_old != remainder_new:
    step += 1
    remainder_new = pow(remainder_old, POWER, MODULO)

    print(f"STEP {step}: {remainder_old}**{POWER} % {MODULO} = {remainder_new}")
    
    remainder_old = remainder_new

step += 1
print(f"STEP {step}: {remainder_old}**{POWER} % {MODULO} = {remainder_new}")

step += 1
remainder_new = pow(remainder_old, POWER, MODULO)
print(f"STEP {step}: {remainder_old}**{POWER} % {MODULO} = {remainder_new}")
print("and so on.")
