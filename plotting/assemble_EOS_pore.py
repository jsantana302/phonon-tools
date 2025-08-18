#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import os
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

def assemble_eos_pore_panel():
    images = ['EOS.png', 'pore_data.png']
    labels = ['a)', 'b)']
    fig, axes = plt.subplots(1, 2, figsize=(16, 5))
    for ax, img, lbl in zip(axes, images, labels):
        if not os.path.exists(img):
            raise FileNotFoundError(f'{img} not found')
        ax.imshow(mpimg.imread(img))
        ax.axis('off')
        ax.text(
            -0.04, 1.02, lbl,
            transform=ax.transAxes,
            ha='left', va='top',
            fontsize=19
        )
    plt.subplots_adjust(
        left=0.05, right=0.95,
        top=0.95, bottom=0.05,
        wspace=0.1
    )
    out_png = 'EOS_pore_panel.png'
    fig.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'✅ Panel saved as {out_png}')

if __name__ == '__main__':
    assemble_eos_pore_panel()

