# Capstan Spring-Brake Model

Python model for a torsional spring-brake assembly using the capstan (belt-friction) effect. Computes static and dynamic performance metrics for a cylindrical spring with friction-based braking — relevant to rotational damping and actuator design.

## What It Models

Given geometry and material parameters, the model computes:

- **Mass and moment of inertia** from cylindrical geometry and material density
- **Torsional spring constant** from material properties
- **Capstan ratio**: Force amplification from friction and wrap angle (F_out / F_in = e^(mu * theta))
- **Brake torque**: Friction torque accounting for wrap angle and normal force
- **Dynamic response**: Time to stop, time to restore, switching frequency
- **Reflected inertia** and net force calculations

## Usage

```bash
python "spring model.py"
```

Default parameters: height = 0.1 m, aluminum density, mu = 0.2.

## Parameters

| Parameter        | Symbol | Default |
|------------------|--------|---------|
| Height           | h      | 0.1 m   |
| Material density | rho    | Al      |
| Friction coeff.  | mu     | 0.2     |
| Wrap angle       | theta  | varies  |

Output prints computed torque, forces, capstan ratio, and brake performance metrics.

## License

MIT
