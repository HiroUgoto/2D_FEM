import matplotlib.pyplot as plt

#--------------------------------------------------------#
def plot_mesh(fem,amp=1.0):
    fig,ax = plt.subplots()

    x = [node.xyz[0] for node in fem.nodes]
    z = [node.xyz[1] for node in fem.nodes]

    ax.set_xlim([min(x)-10,max(x)+10])
    ax.set_ylim([max(z)+5,min(z)-5])
    ax.set_aspect('equal')

    for element in fem.elements:
        f0 = (element.nodes[0].xyz[0], element.nodes[0].xyz[1])
        f1 = (element.nodes[1].xyz[0], element.nodes[1].xyz[1])
        f2 = (element.nodes[2].xyz[0], element.nodes[2].xyz[1])
        f3 = (element.nodes[3].xyz[0], element.nodes[3].xyz[1])
        fpoly = plt.Polygon((f0,f1,f2,f3),ec="k",fc="gray")
        ax.add_patch(fpoly)

    rc = 0.2
    for node in fem.nodes:
        p = plt.Circle((node.xyz[0],node.xyz[1]),rc,color="k")
        ax.add_patch(p)

    plt.show()

#--------------------------------------------------------#
def plot_mesh_update_init():
    _,ax = plt.subplots(figsize=(12,4))
    ax.set_axisbelow(True)
    return ax

def plot_mesh_update(ax,fem,amp=1.0):
    ax.cla()
    ax.grid()

    x = [node.xyz[0] for node in fem.nodes]
    z = [node.xyz[1] for node in fem.nodes]

    ax.set_xlim([min(x)-10,max(x)+10])
    ax.set_ylim([max(z)+5,min(z)-5])
    ax.set_aspect('equal')

    for element in fem.elements:
        f0 = (element.nodes[0].xyz[0]+element.nodes[0].u[0]*amp, element.nodes[0].xyz[1]+element.nodes[0].u[1]*amp)
        f1 = (element.nodes[1].xyz[0]+element.nodes[1].u[0]*amp, element.nodes[1].xyz[1]+element.nodes[1].u[1]*amp)
        f2 = (element.nodes[2].xyz[0]+element.nodes[2].u[0]*amp, element.nodes[2].xyz[1]+element.nodes[2].u[1]*amp)
        f3 = (element.nodes[3].xyz[0]+element.nodes[3].u[0]*amp, element.nodes[3].xyz[1]+element.nodes[3].u[1]*amp)

        fpoly = plt.Polygon((f0,f1,f2,f3),ec="k",fc="gray")
        ax.add_patch(fpoly)

    rc = 0.2
    for node in fem.nodes:
        p = plt.Circle((node.xyz[0]+node.u[0]*amp,node.xyz[1]+node.u[1]*amp),rc,color="k")
        ax.add_patch(p)

    plt.pause(0.001)
