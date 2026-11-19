// File-based testbench for numeriseur_top.
//   input.hex : one line per clock, 48 bits = {x3,x2,x1,x0}, each 12-bit two's complement
//   +tt=<hex16> +win=<dec>   : truth table and coincidence window, loaded through the register interface
// Log (sim/rtl_log.txt):  "C ch n yp yc"  CFD hit of channel ch (n = sample index, y[n-1], y[n])
//                         "H i mask"      hit mask visible after clock i
//                         "T i"           trigger pulse visible after clock i
//                         "S a b c d trig" final scaler values
`timescale 1ns/1ps
module tb_numeriseur;
    parameter integer NS   = 40000;
    parameter integer D    = 3;
    parameter integer K    = 5;
    parameter integer THR  = 90;
    parameter integer HOLD = 16;
    localparam integer EXTRA = 16;

    reg clk = 0;
    reg rst = 1;
    reg signed [11:0] x0 = 0, x1 = 0, x2 = 0, x3 = 0;
    reg tt_we = 0, win_we = 0;
    reg [15:0] tt_data = 0;
    reg [7:0]  win_data = 0;
    wire [3:0] hit;
    wire trig;
    wire [31:0] n0, n1, n2, n3, cnt_a, cnt_b, cnt_c, cnt_d, cnt_trig;
    wire signed [15:0] yp0, yp1, yp2, yp3, yc0, yc1, yc2, yc3;
    reg [47:0] mem [0:NS-1];
    reg [15:0] tt_v;
    integer win_v, i, fd;

    numeriseur_top #(.W(12), .D(D), .K(K), .THR(THR), .HOLD(HOLD), .WW(8)) dut (
        .clk(clk), .rst(rst), .x0(x0), .x1(x1), .x2(x2), .x3(x3),
        .tt_we(tt_we), .tt_data(tt_data), .win_we(win_we), .win_data(win_data),
        .hit(hit), .trig(trig), .n0(n0), .n1(n1), .n2(n2), .n3(n3),
        .yp0(yp0), .yp1(yp1), .yp2(yp2), .yp3(yp3), .yc0(yc0), .yc1(yc1), .yc2(yc2), .yc3(yc3),
        .cnt_a(cnt_a), .cnt_b(cnt_b), .cnt_c(cnt_c), .cnt_d(cnt_d), .cnt_trig(cnt_trig));

    always #2 clk = ~clk;                          // 250 MHz

    initial begin
        if (!$value$plusargs("tt=%h", tt_v)) tt_v = 16'h8888;
        if (!$value$plusargs("win=%d", win_v)) win_v = 3;
        $readmemh("input.hex", mem);
        fd = $fopen("rtl_log.txt", "w");
        repeat (4) @(negedge clk);
        tt_data = tt_v; tt_we = 1; win_data = win_v; win_we = 1;     // load configuration while in reset
        @(negedge clk);
        tt_we = 0; win_we = 0;
        rst = 0;
        for (i = 0; i < NS + EXTRA; i = i + 1) begin
            if (i < NS) {x3, x2, x1, x0} = mem[i]; else {x3, x2, x1, x0} = 48'd0;
            @(negedge clk);
            if (hit[0]) $fwrite(fd, "C 0 %0d %0d %0d\n", n0, yp0, yc0);
            if (hit[1]) $fwrite(fd, "C 1 %0d %0d %0d\n", n1, yp1, yc1);
            if (hit[2]) $fwrite(fd, "C 2 %0d %0d %0d\n", n2, yp2, yc2);
            if (hit[3]) $fwrite(fd, "C 3 %0d %0d %0d\n", n3, yp3, yc3);
            if (hit != 0) $fwrite(fd, "H %0d %0d\n", i, hit);
            if (trig)     $fwrite(fd, "T %0d\n", i);
        end
        $fwrite(fd, "S %0d %0d %0d %0d %0d\n", cnt_a, cnt_b, cnt_c, cnt_d, cnt_trig);
        $fclose(fd);
        $finish;
    end
endmodule
