import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import interact, FloatSlider, Dropdown


# FUNCTION TO CALCULATE FIN TEMPERATURE

def fin_analysis(k=200, h=25, P=0.04, A_c=0.0004, L=0.1, T_b=100, T_inf=25,
                 fin_type="rectangular", tip_condition="insulated"):

    m = np.sqrt(h * P / (k * A_c))
    x = np.linspace(0, L, 200)
    
    # Temperature distribution
    if fin_type == "rectangular":
        if tip_condition == "insulated":
            T_x = T_inf + (T_b - T_inf) * np.cosh(m * (L - x)) / np.cosh(m * L)
            q_fin = np.sqrt(h * P * k * A_c) * (T_b - T_inf) * np.tanh(m * L)
        elif tip_condition == "convective":
            Bi = h * A_c / (k * P)
            T_x = T_inf + (T_b - T_inf) * (np.cosh(m * (L - x)) + Bi/m * np.sinh(m * (L - x))) \
                  / (np.cosh(m * L) + Bi/m * np.sinh(m * L))
            q_fin = k * A_c * m * (T_b - T_inf) * (np.sinh(m * L) + Bi/m * np.cosh(m * L)) \
                    / (np.cosh(m * L) + Bi/m * np.sinh(m * L))
        elif tip_condition == "infinite":
            T_x = T_inf + (T_b - T_inf) * np.exp(-m * x)
            q_fin = np.sqrt(h * P * k * A_c) * (T_b - T_inf)
    elif fin_type == "triangular":
        T_x = T_inf + (T_b - T_inf) * (1 - x/L)
        q_fin = h * P * (T_b - T_inf) * L / 2
    elif fin_type == "pin":
        r = np.sqrt(A_c/np.pi)  # radius of cylindrical pin
        T_x = T_inf + (T_b - T_inf) * np.exp(-np.sqrt(2*h/(k*r)) * x)
        q_fin = 2*np.pi*r*k*(T_b - T_inf) * (1 - np.exp(-np.sqrt(2*h/(k*r)) * L))
    
    # Fin efficiency and effectiveness
    eta_f = q_fin / (h * P * L * (T_b - T_inf))
    epsilon_f = q_fin / (h * P * L * (T_b - T_inf))  # approx. for non-uniform fins

    # Plot
    plt.figure(figsize=(8,5))
    plt.plot(x, T_x, 'r-', linewidth=2)
    plt.xlabel("Length along fin (m)")
    plt.ylabel("Temperature (°C)")
    plt.title(f"{fin_type.capitalize()} Fin ({tip_condition.capitalize()} Tip) Temperature Distribution")
    plt.grid(True)
    plt.show()

    # Print results
    print(f"Heat transfer from fin: {q_fin:.2f} W")
    print(f"Fin efficiency: {eta_f:.3f}")
    print(f"Fin effectiveness: {epsilon_f:.3f}")


# INTERACTIVE SLIDERS AND DROPDOWNS

interact(fin_analysis,
         k=FloatSlider(min=10, max=400, step=5, value=200, description='k (W/m·K)'),
         h=FloatSlider(min=5, max=200, step=5, value=25, description='h (W/m²·K)'),
         P=FloatSlider(min=0.01, max=0.1, step=0.005, value=0.04, description='Perimeter (m)'),
         A_c=FloatSlider(min=0.0001, max=0.01, step=0.0001, value=0.0004, description='Area (m²)'),
         L=FloatSlider(min=0.01, max=0.5, step=0.01, value=0.1, description='Length (m)'),
         T_b=FloatSlider(min=50, max=200, step=5, value=100, description='T_base (°C)'),
         T_inf=FloatSlider(min=0, max=50, step=1, value=25, description='T_ambient (°C)'),
         fin_type=Dropdown(options=["rectangular", "triangular", "pin"], value="rectangular", description='Fin Type'),
         tip_condition=Dropdown(options=["insulated", "convective", "infinite"], value="insulated", description='Tip Condition')
        )
