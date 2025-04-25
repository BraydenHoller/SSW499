import cv2
import os

# Path to the input video file (change this to your file)
video_path = r'C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Plots & Pics & Vids\thickness_permafrost_animation_SSW.mp4'

# Directory where frames will be saved
output_dir = r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Plots & Pics & Vids\perma_frames_SSW"
os.makedirs(output_dir, exist_ok=True)

# Open the video file using OpenCV
cap = cv2.VideoCapture(video_path)
if not cap.isOpened():
    raise IOError(f"Cannot open video file: {video_path}")

frame_count = 0

while True:
    # Read the next frame from the video
    ret, frame = cap.read()
    if not ret:
        break  # End of video

    # Construct a filename for the frame
    frame_filename = os.path.join(output_dir, f'frame_{frame_count:05d}.png')
    
    # Save the frame as an image file (you can choose .jpg, .png, etc.)
    cv2.imwrite(frame_filename, frame)
    
    frame_count += 1

cap.release()
print(f"Extracted {frame_count} frames to the '{output_dir}' directory.")
