using Quon

rz = tait_rz(im*pi)
rx = tait_rx(im*pi)
plot(rz)
plot(rx)
contract!(rz, rx)
plot(rz)