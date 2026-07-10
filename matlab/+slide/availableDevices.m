function devices = availableDevices()
%AVAILABLEDEVICES Return compiled backends that are usable in this runtime.
devices = slide_mex('devices');
devices = reshape(string(devices), 1, []);
end
