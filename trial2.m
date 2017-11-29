% assume the window size is 2 seconds
% there are overlaps in windowing
% assume the step of shif is 1 second

[wav, fs] = audioread('vowels.wav');
wav = wav / max(max(wav));
window_length = 2 * fs;
step = 1 * fs; % has overlap
frame_num = floor((length(wav)-window_length)/step) + 1;

energy = zeros(frame_num, 1);
pos = 1;

for i=1:frame_num
   wav_window = wav(pos:pos + window_length-1);
   energy(i) = 1/window_length * sum(wav_window.^2);
   pos = pos + step;
end
plot(energy);