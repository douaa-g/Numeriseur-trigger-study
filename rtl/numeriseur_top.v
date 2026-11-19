// Four-channel digitiser-trigger core: 4 x digital CFD (rtl/cfd.v) -> trigger4 (coincidence + truth table + scalers).
// One sample per channel per clock (the 250 MS/s ADC clock domain).  The truth table and the window are written
// through a minimal register interface (tt_we / win_we) and are NOT cleared by rst, so they can be loaded while
// the processing is held in reset.  This is a behavioural study block, not a board design.
module numeriseur_top #(
    parameter integer W    = 12,
    parameter integer D    = 3,      // CFD delay in samples (12 ns at 250 MS/s)
    parameter integer K    = 5,      // CFD fraction K/8
    parameter integer THR  = 90,     // arming threshold, LSB
    parameter integer HOLD = 16,
    parameter integer WW   = 8
) (
    input  wire                 clk,
    input  wire                 rst,
    input  wire signed [W-1:0]  x0, x1, x2, x3,
    input  wire                 tt_we,
    input  wire [15:0]          tt_data,
    input  wire                 win_we,
    input  wire [WW-1:0]        win_data,
    output wire [3:0]           hit,
    output wire                 trig,
    output wire [31:0]          n0, n1, n2, n3,
    output wire signed [W+3:0]  yp0, yp1, yp2, yp3, yc0, yc1, yc2, yc3,
    output wire [31:0]          cnt_a, cnt_b, cnt_c, cnt_d, cnt_trig
);
    reg [15:0]   tt  = 16'h0000;
    reg [WW-1:0] win = 0;
    always @(posedge clk) begin
        if (tt_we)  tt  <= tt_data;
        if (win_we) win <= win_data;
    end

    cfd #(.W(W), .D(D), .K(K), .THR(THR), .HOLD(HOLD)) u0 (.clk(clk), .rst(rst), .x(x0), .valid(hit[0]), .n_out(n0), .yp_out(yp0), .yc_out(yc0));
    cfd #(.W(W), .D(D), .K(K), .THR(THR), .HOLD(HOLD)) u1 (.clk(clk), .rst(rst), .x(x1), .valid(hit[1]), .n_out(n1), .yp_out(yp1), .yc_out(yc1));
    cfd #(.W(W), .D(D), .K(K), .THR(THR), .HOLD(HOLD)) u2 (.clk(clk), .rst(rst), .x(x2), .valid(hit[2]), .n_out(n2), .yp_out(yp2), .yc_out(yc2));
    cfd #(.W(W), .D(D), .K(K), .THR(THR), .HOLD(HOLD)) u3 (.clk(clk), .rst(rst), .x(x3), .valid(hit[3]), .n_out(n3), .yp_out(yp3), .yc_out(yc3));

    trigger4 #(.WW(WW)) u_trg (.clk(clk), .rst(rst), .hit(hit), .tt(tt), .win(win), .trig(trig),
                               .cnt_a(cnt_a), .cnt_b(cnt_b), .cnt_c(cnt_c), .cnt_d(cnt_d), .cnt_trig(cnt_trig));
endmodule
