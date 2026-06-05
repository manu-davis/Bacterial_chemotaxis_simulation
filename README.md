A 2D computational biophysics simulation designed to model chemotaxis and physical collision of rod shaped bacteria
---
#chemotaxis

Particles executes a run and tumble motion based on local gradient change

when swimming towards higher concentrations, tumbling reduces and speed increases

when swimming towards lower concentrations, tumbling increases and velocity decreases

#collision(contains a target bacteria and attacking bacteria)

Filters out collision pairs using "scipy.spatial.distance.cdist" and removes considering self collision using an inverted identity matrix mask

Transforms local frame coordinates with target bacteria's centre as origin and rotating such that target lies on x-axis

Collision at centre of target causes it to move by a distance
$$\text{Distance moved} = 2 \cdot dt \cdot \left( \frac{v_{\text{attacker}} - 10}{30} \right) \cdot \cos(\Delta\theta)$$

Collision at any part other than centre causes the target to rotate by an angle
$$\Delta\theta_{\text{target}} = 2 \cdot dt \cdot (-1 \cdot {\text{distance between collision point and centre of target}}) \cdot \left( \frac{v_{\text{attacker}} - 10}{30} \right) \cdot \sin(\Delta\theta)$$

Collision causes attacker bacteria to slow down by a factor
$$\text{factor} = 0.5 \cdot \sin(\Delta\theta)$$

Collision causes attacker to change its dirction slightly towards direction of target
$$\Delta\theta_{\text{collision}} = \Big( ( \Delta\theta + \pi ) \pmod{2\pi} \Big) - \pi$$
