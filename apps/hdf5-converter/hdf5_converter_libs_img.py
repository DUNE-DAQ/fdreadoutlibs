import numpy as np
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import time
import os
import argparse
import warnings
import gc


def create_channel_map_array(drift_direction=0):
    '''
    :param drift_direction: 0 for horizontal drift, 1 for vertical drift
    :return: channel map array
    '''
    # create the channel map array
    if drift_direction == 0:
        channel_map = np.empty((2560, 2), dtype=int)
        channel_map[:, 0] = np.linspace(0, 2559, 2560, dtype=int)
        channel_map[:, 1] = np.concatenate((np.zeros(800), np.ones(800), np.ones(960)*2))        
    elif drift_direction == 1:
        channel_map = np.empty((3072, 2), dtype=int)
        channel_map[:, 0] = np.linspace(0, 3071, 3072, dtype=int)
        channel_map[:, 1] = np.concatenate((np.zeros(952), np.ones(952), np.ones(1168)*2))
    elif drift_direction == 2:
        channel_map = np.empty((128, 2), dtype=int)
        channel_map[:, 0] = np.linspace(0, 127, 128, dtype=int)
        channel_map[:, 1] = np.concatenate((np.zeros(39), np.ones(40), np.ones(49)*2))

    return channel_map


def group_maker(all_tps, channel_map, ticks_limit=100, channel_limit=20, min_tps_to_group=2):
    '''
    :param all_tps: all trigger primitives in the event
    :param channel_map: channel map
    :param ticks_limit: maximum time window to consider
    :param channel_limit: maximum number of channels to consider
    :param min_tps_to_group: minimum number of hits to consider
    :return: list of groups
    '''
    # if i have tps from multiple planes, i group only by time
    total_channels = channel_map.shape[0]
    if np.unique(channel_map[all_tps[:, 3]% total_channels, 1]).shape[0] != 1:
        print('Warning: multiple planes in the event. Grouping only by time.')
        return group_maker_only_by_time(all_tps, channel_map, ticks_limit=ticks_limit, channel_limit=channel_limit, min_tps_to_group=min_tps_to_group)


    # create a new algorithm to group the TPs
    groups = []
    buffer = []
    # loop over the TPs
    for tp in all_tps:
        if len(buffer)==0:
            buffer.append([tp])
        else:
            buffer_copy = buffer.copy()
            buffer = []
            appended = False
            idx = 0
            for candidate in (buffer_copy):
                time_cond = (tp[0] - np.max(np.array(candidate)[:, 0]+np.array(candidate)[:,1])) < ticks_limit
                if time_cond:
                    chan_cond = np.min(np.abs(tp[3] - np.array(candidate)[:, 3])) < channel_limit
                    if chan_cond:
                        if not appended:
                            candidate.append(tp)
                            buffer.append(candidate)
                            idx_appended = idx
                            appended = True
                        else:
                            for tp2 in candidate:
                                buffer[idx_appended].append(tp2)
                    else:
                        buffer.append(candidate)
                        idx += 1
                else:
                    if len(candidate) >= min_tps_to_group:
                        groups.append(np.array(candidate, dtype=int))
            if not appended:
                buffer.append([tp])
    if len(buffer) > 0:
        for candidate in buffer:
            if len(candidate) >= min_tps_to_group:
                groups.append(np.array(candidate, dtype=int))

    groups = np.array(groups, dtype=object)
    return groups



def group_maker_only_by_time(all_tps, channel_map, ticks_limit=100, channel_limit=20, min_tps_to_group=4):
    n_numbers = all_tps.shape[1]
    groups = []
    current_group = np.array([all_tps[i] for i in range(n_numbers)])
    for tp in all_tps[1:]:
        tp = [tp[i] for i in range(n_numbers)]
        if (tp[0] - np.max(current_group[:, 0]+current_group[:,1])) < ticks_limit:
            current_group = np.append(current_group, [tp], axis=0)
        else:
            if len(current_group) >= min_tps_to_group:
                groups.append(current_group)
            current_group = np.array([tp])
            

    return groups

def from_tp_to_imgs(tps, make_fixed_size=False, width=500, height=1000, x_margin=10, y_margin=100, y_min_overall=-1, y_max_overall=-1):
    '''
    :param tps: all trigger primitives to draw
    :param make_fixed_size: if True, the image will have fixed size, otherwise it will be as big as the TPs
    :param width: width of the image
    :param height: height of the image
    :param x_margin: margin on the x axis
    :param y_margin: margin on the y axis
    :param y_min_overall: minimum y value of the image
    :param y_max_overall: maximum y value of the image
    :return: image
    '''

    if y_min_overall != -1:
        t_start = y_min_overall
    else:
        t_start = tps[0, 0] 

    if y_max_overall != -1:
        t_end = y_max_overall
    else:
        t_end = np.max(tps[:, 0] + tps[:, 1])
    
    x_max = (tps[:, 3].max())
    x_min = (tps[:, 3].min())

    x_range = x_max - x_min 
    y_range = int((t_end - t_start))

    # create the image
    if make_fixed_size:
        img_width = width
        img_height = height
    else:
        img_width =  x_range + 2*x_margin
        img_height =  y_range + 2*y_margin
    img = np.zeros((img_height, img_width))
    # fill image
    if (not make_fixed_size):
        for tp in tps:
            x = (tp[3] - x_min) + x_margin
            y_start = (tp[0] - t_start) + y_margin
            y_end = (tp[0] + tp[1] - t_start) + y_margin
            img[int(y_start)-1:int(y_end), int(x)-1] = tp[4]/(y_end - y_start)


    else:
    # We stretch the image inwards if needed but we do not upscale it. In this second case we build a padded image
        if img_width < x_range:
            stretch_x = True
            print('Warning: image width is smaller than the range of the TPs. The image will be stretched inwards.')

        else:
            x_margin = (img_width - x_range)/2
            stretch_x = False
        
        if img_height < y_range:
            stretch_y = True
            print('Warning: image height is smaller than the range of the TPs. The image will be stretched inwards.')

        else:
            y_margin = (img_height - y_range)/2
            stretch_y = False


        if stretch_x & stretch_y:
            for tp in tps:
                x=(tp[3] - x_min)/x_range * (img_width - 2*x_margin) + x_margin
                y_start = (tp[0] - t_start)/y_range * (img_height - 2*y_margin) + y_margin
                y_end = (tp[0] + tp[1] - t_start)/y_range * (img_height - 2*y_margin) + y_margin
                img[int(y_start)-1:int(y_end), int(x)-1] = tp[4]/(y_end - y_start)
        elif stretch_x:
            for tp in tps:                
                x=(tp[3] - x_min)/x_range * (img_width - 2*x_margin) + x_margin
                y_start = (tp[0] - t_start) + y_margin
                y_end = (tp[0] + tp[1] - t_start) + y_margin
                img[int(y_start)-1:int(y_end), int(x)-1] = tp[4]/(y_end - y_start)
        elif stretch_y:
            for tp in tps:
                x = (tp[3] - x_min) + x_margin
                y_start = (tp[0] - t_start)/y_range * (img_height - 2*y_margin) + y_margin
                y_end = (tp[0] + tp[1] - t_start)/y_range * (img_height - 2*y_margin) + y_margin
                img[int(y_start):int(y_end), int(x)-1] = tp[4]/(y_end - y_start)

        else:
            for tp in tps:
                x = (tp[3] - x_min) + x_margin
                y_start = (tp[0] - t_start) + y_margin
                y_end = (tp[0] + tp[1] - t_start) + y_margin
                img[int(y_start)-1:int(y_end), int(x)-1] = tp[4]/(y_end - y_start)
   
    return img


def all_views_img_maker(tps, channel_map, min_tps_to_create_img=2, make_fixed_size=False, width=500, height=1000, x_margin=10, y_margin=200):
    '''
    :param tps: all trigger primitives to draw
    :param channel_map: channel map array
    :param min_tps_to_create_img: minimum number of TPs to create an image
    :param make_fixed_size: if True, the image will have fixed size, otherwise it will be as big as the TPs
    :param width: width of the image
    :param height: height of the image
    :param x_margin: margin on the x axis
    :param y_margin: margin on the y axis
    :return: images or -1 if there are not enough TPs
    '''

    total_channels = channel_map.shape[0]
    # U plane, take only the tps where the corrisponding position in the channel map is 0      
    tps_u = tps[np.where(channel_map[tps[:, 3]% total_channels, 1] == 0)]

    # V plane, take only the tps where the corrisponding position in the channel map is 1
    tps_v = tps[np.where(channel_map[tps[:, 3]% total_channels, 1] == 1)]

    # X plane, take only the tps where the corrisponding position in the channel map is 2
    tps_x = tps[np.where(channel_map[tps[:, 3]% total_channels, 1] == 2)]
    print(tps_u.shape)
    print(tps_v.shape)
    print(tps_x.shape)

    img_u, img_v, img_x = np.array([[-1]]), np.array([[-1]]), np.array([[-1]])

    y_min_overall = tps[0,0]
    y_max_overall = np.max(tps[:,0] + tps[:,1])

    if tps_u.shape[0] >= min_tps_to_create_img:
        img_u = from_tp_to_imgs(tps_u, make_fixed_size=make_fixed_size, width=width, height=height, x_margin=x_margin, y_margin=y_margin, y_min_overall=y_min_overall, y_max_overall=y_max_overall)

    if tps_v.shape[0] >= min_tps_to_create_img: 
        img_v = from_tp_to_imgs(tps_v, make_fixed_size=make_fixed_size, width=width, height=height, x_margin=x_margin, y_margin=y_margin, y_min_overall=y_min_overall, y_max_overall=y_max_overall)

    if tps_x.shape[0] >= min_tps_to_create_img:
        img_x = from_tp_to_imgs(tps_x, make_fixed_size=make_fixed_size, width=width, height=height, x_margin=x_margin, y_margin=y_margin, y_min_overall=y_min_overall, y_max_overall=y_max_overall)


    return img_u, img_v, img_x


def save_img(all_TPs, channel_map,save_path, outname='test', min_tps_to_create_img=2, make_fixed_size=False, width=500, height=1000, x_margin=10, y_margin=200):
    '''
    :param all_TPs: all trigger primitives in the event
    :param channel_map: channel map
    :param save_path: path where to save the images
    :param outname: base name for the images
    :param min_tps: minimum number of TPs to create the image
    :param make_fixed_size: if True, the image will have fixed size, otherwise it will be as big as the TPs
    :param width: width of the image if make_fixed_size is True
    :param height: height of the image if make_fixed_size is True
    :param x_margin: margin on the x axis if make_fixed_size is False
    :param y_margin: margin on the y axis if make_fixed_size is False
    '''
    total_channels = channel_map.shape[0]

    #create images
    img_u, img_v, img_x = all_views_img_maker(all_TPs, channel_map, min_tps_to_create_img=min_tps_to_create_img, make_fixed_size=make_fixed_size, width=width, height=height, x_margin=x_margin, y_margin=y_margin)
    max_pixel_value_overall = np.max([np.max(img_u), np.max(img_v), np.max(img_x)])

    #save images
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    x_max = (all_TPs[:, 3].max())
    x_min = (all_TPs[:, 3].min())

    x_range = x_max - x_min 
    t_start = all_TPs[0, 0] 

    t_end = np.max(all_TPs[:, 0] + all_TPs[:, 1])
    y_range = int((t_end - t_start))

    # create the image
    if make_fixed_size:
        img_width = width
        img_height = height
        if img_width > x_range:
            x_margin = (img_width - x_range)/2    
        if img_height > y_range:
            y_margin = (img_height - y_range)/2
    else:
        img_width =  x_range + 2*x_margin
        img_height =  y_range + 2*y_margin
    yticks_labels = [t_start-y_margin + i*(y_range + 2*y_margin)//10 for i in range(10)]

    n_views = 0

    if img_u[0, 0] != -1:
        tps_u = all_TPs[np.where(channel_map[all_TPs[:, 3]% total_channels, 1] == 0)]
        
        x_min_u = (tps_u[:, 3].min())
        x_max_u = (tps_u[:, 3].max())
        x_range_u = x_max_u - x_min_u

        x_margin_u = x_margin

        if make_fixed_size:
            if img_width > x_range_u:
                x_margin_u = (img_width - x_range_u)/2    
        xticks_labels_u = [x_min_u-x_margin_u + i*(x_range_u + 2*x_margin_u)//2 for i in range(2)]

        n_views += 1
        plt.figure(figsize=(10, 26))
        plt.title('U plane')
        plt.imshow(img_u)
        # plt.colorbar()
        # add x and y labels
        plt.xlabel("Channel")
        plt.ylabel("Time (ticks)")
        # set y axis ticks
        plt.yticks(ticks=np.arange(0, img_v.shape[0], img_v.shape[0]/10), labels=yticks_labels)
        # set x axis ticks
        plt.xticks(ticks=np.arange(0, img_u.shape[1], img_u.shape[1]/2), labels=xticks_labels_u)

        # save the image, with a bbox in inches smaller than the default but bigger than tight
        plt.savefig(save_path+ 'u_' + os.path.basename(outname) + '.png', bbox_inches='tight', pad_inches=1)
        plt.close()

    if img_v[0, 0] != -1:
        tps_v = all_TPs[np.where(channel_map[all_TPs[:, 3]% total_channels, 1] == 1)]

        x_min_v = (tps_v[:, 3].min())
        x_max_v = (tps_v[:, 3].max())
        x_range_v = x_max_v - x_min_v

        x_margin_v = x_margin

        if make_fixed_size:
            if img_width > x_range_v:
                x_margin_v = (img_width - x_range_v)/2
        xticks_labels_v = [x_min_v-x_margin_v + i*(x_range_v + 2*x_margin_v)//2 for i in range(2)]
        
        n_views += 1
        plt.figure(figsize=(10, 26))    
        plt.title('V plane')
        plt.imshow(img_v)
        # plt.colorbar()
        # add x and y labels
        plt.xlabel("Channel")
        plt.ylabel("Time (ticks)")

        # set y axis ticks
        plt.yticks(ticks=np.arange(0, img_v.shape[0], img_v.shape[0]/10), labels=yticks_labels)
        # set x axis ticks
        plt.xticks(ticks=np.arange(0, img_v.shape[1], img_v.shape[1]/2), labels=xticks_labels_v)

        # save the image, with a bbox in inches smaller than the default but bigger than tight
        plt.savefig(save_path+ 'v_' + os.path.basename(outname) + '.png', bbox_inches='tight', pad_inches=1)
        plt.close()

    if img_x[0, 0] != -1:
        tps_x = all_TPs[np.where(channel_map[all_TPs[:, 3]% total_channels, 1] == 2)]

        x_min_x = (tps_x[:, 3].min())
        x_max_x = (tps_x[:, 3].max())
        x_range_x = x_max_x - x_min_x

        x_margin_x = x_margin

        if make_fixed_size:
            if img_width > x_range_x:
                x_margin_x = (img_width - x_range_x)/2
        xticks_labels_x = [x_min_x-x_margin_x + i*(x_range_x + 2*x_margin_x)//2 for i in range(2)]

        n_views += 1
        plt.figure(figsize=(10, 26))
        plt.title('X plane')
        plt.imshow(img_x)
        # plt.colorbar()
        # add x and y labels
        plt.xlabel("Channel")
        plt.ylabel("Time (ticks)")
        # set y axis ticks
        plt.yticks(ticks=np.arange(0, img_v.shape[0], img_v.shape[0]/10), labels=yticks_labels)
        # set x axis ticks
        plt.xticks(ticks=np.arange(0, img_x.shape[1], img_x.shape[1]/2), labels=xticks_labels_x)

        # save the image, with a bbox in inches smaller than the default but bigger than tight
        plt.savefig(save_path+ 'x_' + os.path.basename(outname) + '.png', bbox_inches='tight', pad_inches=1)
        plt.close()
    if n_views == 0:
        print(f'No images saved! Do a better grouping algorithm! You had {all_TPs.shape[0]} tps and {min_tps_to_create_img} as min_tps_to_create_img.' )
        print(f'You have {all_TPs[np.where(channel_map[all_TPs[:, 3]% total_channels, 1] == 0)].shape[0]} U tps, {all_TPs[np.where(channel_map[all_TPs[:, 3]% total_channels, 1] == 1)].shape[0]} V tps and {all_TPs[np.where(channel_map[all_TPs[:, 3]% total_channels, 1] == 2)].shape[0]} Z tps.')


    if n_views > 1:
        fig = plt.figure(figsize=(10, 26))
        grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                        nrows_ncols=(1,3),
                        axes_pad=0.5,
                        share_all=True,
                        cbar_location="right",
                        cbar_mode="single",
                        cbar_size="30%",
                        cbar_pad=0.25,
                        )   


        if img_u[0, 0] != -1:
            im = grid[0].imshow(img_u, vmin=0, vmax=max_pixel_value_overall)
            grid[0].set_title('U plane')
            # grid[0].set_xticks(np.arange(0, img_u.shape[1], img_u.shape[1]/2))
            # grid[0].set_xticklabels(xticks_labels_u)
        if img_v[0, 0] != -1:
            im = grid[1].imshow(img_v, vmin=0, vmax=max_pixel_value_overall)
            grid[1].set_title('V plane')
            # grid[1].set_xticks(np.arange(0, img_v.shape[1], img_v.shape[1]/2))
            # grid[1].set_xticklabels(xticks_labels_v)
        if img_x[0, 0] != -1:
            im = grid[2].imshow(img_x, vmin=0, vmax=max_pixel_value_overall)
            grid[2].set_title('X plane')
            # grid[2].set_xticks(np.arange(0, img_x.shape[1], img_x.shape[1]/2))
            # grid[2].set_xticklabels(xticks_labels_x)

        grid.cbar_axes[0].colorbar(im)
        # grid.axes_llc.set_yticks(yticks_labels)
        # use the same yticks_labels for all the images
        grid.axes_llc.set_yticks(np.arange(0, img_v.shape[0], img_v.shape[0]/10))
        grid.axes_llc.set_yticklabels(yticks_labels)
        
        # save the image
        plt.savefig(save_path+ 'multiview_' + os.path.basename(outname) + '.png')
        plt.close()

