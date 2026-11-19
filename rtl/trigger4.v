// Four-input configurable trigger: coincidence gate + truth-table lookup + scalers.
//
// Same architecture idea as the LPSC UCTM/UCTM2 modules (arXiv:1603.06816): the user's logic equation is
// evaluated by software for the 16 input combinations and loaded as a 16-bit truth table, so a new equation
// never requires a new synthesis.
//
//  * the first discriminator hit opens a gate; hits arriving up to WIN clock cycles after the first one are
//    collected (sticky bits)
//  * when the gate closes, the collected pattern {D,C,B,A} addresses the truth table and, if the entry is 1,
//    a one-cycle trigger is produced.  The decision is therefore taken once per gate, which is what makes
//    vetoes (A & B & !D) and XOR behave correctly; the price is a fixed latency of WIN+1 cycles
//  * 32-bit scalers count input hits and triggers
// WIN must be >= 1 (0 is treated as 1).
module trigger4 #(
    parameter integer WW = 8        // width of the window register
) (
    input  wire          clk,
    input  wire          rst,
    input  wire [3:0]    hit,       // one-cycle pulses from the four channel discriminators
    input  wire [15:0]   tt,        // truth table, index = {D,C,B,A}
    input  wire [WW-1:0] win,       // gate length in clock cycles after the first hit
    output reg           trig,
    output reg  [31:0]   cnt_a, cnt_b, cnt_c, cnt_d, cnt_trig
);
    reg          open_g;
    reg [WW-1:0] gc;
    reg [3:0]    sticky;
    wire [3:0]   acc = sticky | hit;

    always @(posedge clk) begin
        if (rst) begin
            open_g <= 0; gc <= 0; sticky <= 0; trig <= 0;
            cnt_a <= 0; cnt_b <= 0; cnt_c <= 0; cnt_d <= 0; cnt_trig <= 0;
        end else begin
            trig <= 0;
            if (!open_g) begin
                if (|hit) begin
                    open_g <= 1;
                    sticky <= hit;
                    gc     <= (win == 0) ? 1 : win;
                end
            end else begin
                sticky <= acc;
                if (gc == 1) begin
                    open_g <= 0;
                    sticky <= 0;
                    trig   <= tt[acc];
                    cnt_trig <= cnt_trig + tt[acc];
                end else begin
                    gc <= gc - 1;
                end
            end
            cnt_a <= cnt_a + hit[0];
            cnt_b <= cnt_b + hit[1];
            cnt_c <= cnt_c + hit[2];
            cnt_d <= cnt_d + hit[3];
        end
    end
endmodule
