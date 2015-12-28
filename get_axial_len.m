function axial_len = get_axial_len(a, b, c, d, f, g)

axial_len = zeros(1, 2);
axial_len(1) = sqrt(2 * (a * f^2 + c * d^2 + g * b^2 - 2 * b * d * f - a * c * g) / ( (b^2 - a * c) * (sqrt( (a - c)^2 + 4 * b^2) - (a + c) ) ) );
axial_len(2) = sqrt(2 * (a * f^2 + c * d^2 + g * b^2 - 2 * b * d * f - a * c * g) / ( (b^2 - a * c) * (-sqrt( (a - c)^2 + 4 * b^2) - (a + c) ) ) );