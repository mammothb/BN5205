def psi(num, der, xi_1, xi_2):
    if der == 0:
        if num == 0:
            return (1 - xi_1) * (1 - xi_2)
        elif num == 1:
            return xi_1 * (1 - xi_2)
        elif num == 2:
            return (1 - xi_1) * xi_2
        elif num == 3:
            return xi_1 * xi_2
        else:
            return 0
    elif der == 1:
        if num == 0:
            return -1 * (1 - xi_2)
        elif num == 1:
            return 1 - xi_2
        elif num == 2:
            return -xi_2
        elif num == 3:
            return xi_2
        else:
            return 0
    elif der == 2:
        if num == 0:
            return -1 * (1 - xi_1)
        elif num == 1:
            return -xi_1
        elif num == 2:
            return 1 - xi_1
        elif num == 3:
            return xi_1
        else:
            return 0
    else:
        return 0
