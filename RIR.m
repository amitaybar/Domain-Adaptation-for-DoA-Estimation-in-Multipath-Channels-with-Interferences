function [ h ] = RIR( s,fs,r,n,c,L,order)


SoundVel    = 342;                   % Sound velocity (m/s)
mtype       = 'omnidirectional';  % Type of microphone
RoomDim     = 3;                    % Room dimension
orientation = 0;            % Microphone orientation (rad)
hp_filter   = 1;              % Enable high-pass filter
h           = rir_generator(SoundVel, fs, r, s, L, c, n, mtype,order,RoomDim,orientation,hp_filter);

end
