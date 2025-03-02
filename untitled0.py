import numpy as np
import matplotlib.pyplot as plt

# Inizializzazione degli array
x = np.zeros(1000)
y = np.zeros(1000)
x[0] = 0.5
y[0] = 0.5

# Calcolo delle sequenze per le due equazioni
for n in range(0, len(x) - 1):
    x[n + 1] = 3.9 * x[n] * (1 - x[n])
    y[n + 1] = 3.9001 * y[n] * (1 - y[n])

k = np.arange(len(x))


# Creazione di due grafici separati
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

# Primo grafico per x[n]
ax1.plot(k, x, color='b', label='x[n] = 3.9 * x[n] * (1 - x[n])')
ax1.set_xlabel('n')
ax1.set_ylabel('Valori di x')
ax1.set_title('Equazione x[n] = 3.9 * x[n] * (1 - x[n])')
ax1.grid(True)
ax1.legend()

# Secondo grafico per y[n]
ax2.plot(k, y, color='r', label='y[n] = 3.9001 * y[n] * (1 - y[n])')
ax2.set_xlabel('n')
ax2.set_ylabel('Valori di y')
ax2.set_title('Equazione y[n] = 3.9001 * y[n] * (1 - y[n])')
ax2.grid(True)
ax2.legend()

# Mostra i grafici
plt.tight_layout()
plt.show()
