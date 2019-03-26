if (c == 'Y') {
  o = 0;
  o += feed_forward(s, netphosk_InsR_1, 17, 8);
  o += feed_forward(s, netphosk_InsR_2, 9, 8);
  o += feed_forward(s, netphosk_InsR_3, 7, 4);
  o /= 3;
}
