# Witch of Agnesi

This setup replicates the "Linear Nonhydrostatic Mountain" problem from Giraldo & Restelli 2008,
JCP. Note that this problem is sensitive to the damping layer extent (Rayleigh damping at zhi;
sponge regions at xlo, xhi) as well as the damping coefficient. We tuned the damping coefficient
within the recommended range from Durran & Klemp 1983.

The solver executable is compiled in Exec/DryRegTests/WitchOfAgnesi and provides the analytical
surface geometry.
