def psi(num, der, xi_1, xi_2):
    return {
        0: {
            0: (1 - xi_1) * (1 - xi_2),
            1: xi_1 * (1 - xi_2),
            2: (1 - xi_1) * xi_2,
            3: xi_1 * xi_2
        },
        1: {
            0: -(1 - xi_2),
            1: 1 - xi_2,
            2: -xi_2,
            3: xi_2
        },
        2: {
            0: -(1 - xi_1),
            1: -xi_1,
            2: 1 - xi_1,
            3: xi_1
        }
    }.get(der, 0).get(num, 0)
