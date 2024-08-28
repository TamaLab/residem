import matplotlib.pyplot as plt
import imageio.v2 as imageio
import cv2
import os
import numpy as np


class images2movies_extra:
    def __init__(self, folder_path=None):
        self.folder_path = folder_path if folder_path is not None else os.getcwd()
        self.image_main = os.path.join(folder_path, "diff_map/diff_map.png")
        self.folder_path_2 = os.path.join(folder_path,"extra")
        self.folder_path_3 = os.path.join(folder_path,"extra_diff_map")
        self.imagesave()
        self.imagesave_2()

    def imagesave(self):
        images = [f'%s/extra_occu_%s.png' %
                  (self.folder_path_2, x) for x in range(1, 101)]
        existing_images = [image for image in images if os.path.exists(image)]
        existing_numbers = [(int(image.split("_")[-1].split(".")[0]),image) for image in existing_images]

        # labels = [f'Occupancy -%s' % x[0] for x in existing_numbers]
        os.makedirs(f'{self.folder_path}/extra_image_occu', exist_ok=True)

        #for i in range(len(images)):# images
        for i in existing_numbers:
            fig, axs = plt.subplots()
            fig = plt.figure() #figsize=(25, 10)
            img = plt.imread(i[1])
            plt.imshow(img)
            # plt.title(labels[i],fontsize=12)
            plt.title(f"Occupancy -%s"%i[0], fontsize=12)
            plt.axis('off')
            plt.savefig(f'{self.folder_path}/extra_image_occu/labeled_image{i[0]}.png')
            plt.close()

        with imageio.get_writer(f'{self.folder_path}/movie_extra_occu.gif', mode='I') as writer:
            #for i in range(len(images)):
            for i in existing_numbers:
                image = imageio.imread(
                    f'{self.folder_path}/extra_image_occu/labeled_image{i[0]}.png')
                writer.append_data(image)

        images_to_movie = [f'{self.folder_path}/extra_image_occu/labeled_image{i[0]}.png'for i in existing_numbers]
           # f'{self.folder_path}/extra_image_occu/labeled_image{i}.png' for i in range(len(images))]



        frame = cv2.imread(images_to_movie[0])
        height, width, layers = frame.shape
        video = cv2.VideoWriter(
            f'{self.folder_path}/movie_extra_image.mp4', cv2.VideoWriter_fourcc(*'mp4v'), 1, (width, height))
        for image_x in images_to_movie:
            video.write(cv2.imread(image_x))

        video.release()

    def imagesave_2(self):
        images = [f'%s/extra_diff_map_occu_%s.png' %
                  (self.folder_path_3, x) for x in range(1, 101)]
        existing_images = [image for image in images if os.path.exists(image)]
        existing_numbers = [(int(image.split("_")[-1].split(".")[0]),image) for image in existing_images]

        # labels = [f'Extrapolated difference map occupancy-%s' %
        #           x for x in range(1, 101)]

        # labels = [f'Extrapolated difference map occupancy-%s' %
        #           x[0] for x in existing_numbers]

        os.makedirs(f'%s/extra_diff_image_occu'%self.folder_path, exist_ok=True)
        # for i in range(len(images)):
        for i in existing_numbers:
            fig, axs = plt.subplots(1, 2, constrained_layout=True) #figsize=(8, 10)
            # fig = plt.figure(figsize=(25, 10))
            image_de = plt.imread(self.image_main)
            axs[0].imshow(image_de)
            axs[0].axis('off')
            axs[0].set_title("Isomorphous difference map",fontsize=10)
            # img = plt.imread(images[i])
            img = plt.imread(i[1])
            axs[1].imshow(img)
            axs[1].axis('off')
            # axs[1].set_title(labels[i],fontsize=10)
            axs[1].set_title(f'Extrapolated difference map occupancy-%s'%i[0], fontsize=10)
            plt.tight_layout()
            # plt.subplots_adjust(top=0.95, bottom=0.95)
            plt.savefig(f'{self.folder_path}/extra_diff_image_occu/labeled_image{i[0]}.png')
            plt.close()

        with imageio.get_writer(f'{self.folder_path}/movie_extra_diff_occu.gif', mode='I') as writer:
            # for i in range(len(images)):
            for i in existing_numbers:
                image = imageio.imread(
                    f'{self.folder_path}/extra_diff_image_occu/labeled_image{i[0]}.png')
                writer.append_data(image)

        images_to_movie = [
            f'{self.folder_path}/extra_diff_image_occu/labeled_image{i[0]}.png' for i in existing_numbers]

        frame = cv2.imread(images_to_movie[0])
        height, width, layers = frame.shape
        video = cv2.VideoWriter(
            f'{self.folder_path}/movie_extra_diff_image.mp4', cv2.VideoWriter_fourcc(*'mp4v'), 1, (width, height))
        for image_x in images_to_movie:
            video.write(cv2.imread(image_x))

        video.release()





# frame = cv2.imread(images[0])
# height, width, layers = frame.shape
#
# # Create a VideoWriter object
# video = cv2.VideoWriter(
#     'movie.mp4', cv2.VideoWriter_fourcc(*'mp4v'), 1, (width, height))
#
# # Add images to the video writer
# for image in images:
#     video.write(cv2.imread(image))
#
# # Close the video writer
# video.release()


def run(folder):
    images2movies_extra(folder_path=folder)


def run_all():
    import argparse
    from datetime import datetime

    start = datetime.now()

    parser = argparse.ArgumentParser(prog='residem_image_2_movie', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-f", "--folder", type=str, help="folder path")


    args = parser.parse_args()


    run(args.folder)
    print(datetime.now() - start)


if __name__ == "__main__":
    run_all()


