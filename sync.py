import signal
import sys
import time

import serial

import ika.ik
from autocal.lstsq import calibrate_angles, calibrate_lengths, get_point


def interrupt_handler(signal, frame):
    print(" [interrupt] ")
    sys.exit(0)


if __name__ == "__main__":
    signal.signal(signal.SIGINT, interrupt_handler)

    calibration_points = [
        (0.0, 0.0, 0.0),
        (0.0, 100.0, 0.0),
        (100.0, 100.0, 0.0),
        (100.0, 0.0, 0.0),
    ]

    ser = serial.Serial("/dev/ttyACM0", 115200, timeout=1)
    print("using port", ser.name, "@", ser.baudrate)
    ser2 = serial.Serial("/dev/ttyACM1", 115200, timeout=1)
    print("using port", ser2.name, "@", ser.baudrate)

    coeffs = calibrate_angles(ser)
    res = calibrate_lengths(ser, coeffs, calibration_points)
    if res is not None:
        best_params, base_angle = res
    else:
        quit()

    while True:
        target = get_point(ser, coeffs, best_params, base_angle)
        print(target)
        try:
            j1, j2, j3, j4, j5 = ika.ik.get_joints_to(target[0], target[1], target[2])
        except:
            continue
        ser2.write(f"{j1},{j2},{j3},{j4},{j5}%".encode())
        # ser2.flush()
        time.sleep(0.01)
        ser2.reset_input_buffer()
