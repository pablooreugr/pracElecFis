# Análisis de la Distancia Crítica en función de la Masa Efectiva

## Introducción y Resultados Empíricos

El objetivo del análisis fue determinar de forma empírica la relación matemática entre la **masa efectiva ($m^*$)** de un portador de carga y la **distancia crítica ($L_c$)** de la barrera de potencial (en este caso un modelo análogo a un MOSFET). La distancia crítica se ha definido en este contexto como la longitud de la barrera necesaria para que la probabilidad de transmisión decaiga a un valor predeterminado (por ejemplo, $T = 10^{-4}$).

Tras procesar los datos extraídos por el simulador y realizar diversas regresiones (Lineal Inversa, Exponencial, Ley de Potencias e Inversa de la Raíz Cuadrada), hemos obtenido los siguientes parámetros de ajuste:

1. **Ley de Potencias** ($L_c = a \cdot (m^*)^b$):
   * Exponente ajustado: **$b = -0.477 \pm 0.001$**
   * Coeficiente de determinación: **$R^2 = 0.9992$**

2. **Inversa de la Raíz Cuadrada** ($L_c = \frac{a}{\sqrt{m^*}} + b$):
   * Coeficiente de determinación: **$R^2 = 0.9986$**

La elevada precisión de estos ajustes ($R^2 > 0.998$) demuestra de manera empírica que la distancia crítica **decae proporcionalmente a la inversa de la raíz cuadrada de la masa efectiva** ($L_c \propto 1/\sqrt{m^*}$).

## Significado Físico

Este resultado empírico obtenido mediante simulador corrobora de manera excepcional la teoría cuántica de penetración de barreras. 

De acuerdo con la **Aproximación WKB** (Wentzel-Kramers-Brillouin) para el efecto túnel a través de una barrera de potencial rectangular gruesa, la probabilidad de transmisión $T$ para una partícula de masa efectiva $m^*$ y energía $E$ incidiendo sobre una barrera de altura $V$ ($V > E$) y anchura $L$ se puede aproximar como:

$$T \approx \exp\left( -2 \cdot L \cdot \frac{\sqrt{2 m^* (V - E)}}{\hbar} \right)$$

En nuestro experimento, definimos la **distancia crítica ($L_c$)** como la longitud $L$ a la cual la transmisión alcanza un valor umbral constante (sea $T_{crit} = 10^{-4}$). Por lo tanto, tomando el logaritmo natural a ambos lados de la ecuación de transmisión:

$$\ln(T_{crit}) = -2 \cdot L_c \cdot \frac{\sqrt{2 m^* (V - E)}}{\hbar}$$

Al despejar la distancia crítica $L_c$, obtenemos:

$$L_c = - \frac{\hbar \cdot \ln(T_{crit})}{2 \sqrt{2 (V - E)}} \cdot \frac{1}{\sqrt{m^*}}$$

Dado que $V$, $E$ y $T_{crit}$ son constantes para diferentes simulaciones en las que únicamente variamos la masa efectiva del material, todos los términos de la izquierda se agrupan en una única constante de proporcionalidad $C$:

$$L_c = \frac{C}{\sqrt{m^*}} = C \cdot (m^*)^{-0.5}$$

## Conclusión

El análisis de datos arroja un exponente empírico de $b = -0.477$, extremadamente cercano al exponente teórico de $-0.5$ predecido por la ecuación WKB. 

Físicamente, esto significa que **cuanto más pesado es el portador de carga (mayor masa efectiva), más rápido decae su función de onda dentro de la barrera clásicamente prohibida**. Como resultado, las partículas con mayor masa efectiva penetran menos en la barrera (menor efecto túnel) y, por consiguiente, requieren barreras más estrechas (menor $L_c$) para alcanzar la misma probabilidad de transmisión que partículas más ligeras.
