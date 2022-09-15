N_TEETH = int(input())
X_LAST_TOOTH = (N_TEETH - 1) * 3 + 5

print(f'0.4\n1\n{4 * N_TEETH + 5}\n4 3 4 5\n4 5', end='') # learning rate, pull attraction, number of segments

for i in range(5, X_LAST_TOOTH, 3):
    print(f' {i} 9\n{i} 9 {i + 1} 5\n{i + 1} 5 {i + 2} 5\n{i + 2} 5', end='')

print(f' {X_LAST_TOOTH} 9\n{X_LAST_TOOTH} 9 {X_LAST_TOOTH + 1} 5\n{X_LAST_TOOTH + 1} 5 {X_LAST_TOOTH + 1} 3\n{X_LAST_TOOTH + 1} 3 4 3\n{N_TEETH}')

for i in range(1, N_TEETH + 1):
    print(f'{4 + 0.1 * i} 4')

