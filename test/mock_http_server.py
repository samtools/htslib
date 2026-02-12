#!/usr/bin/env python3
"""Mock HTTP server for testing hfile_libcurl retry logic.

Usage:
    python test/mock_http_server.py --mode <mode> --file <path> [--fail-count N] [--port 0]

Prints the allocated port number to stdout on startup, then serves requests.

Modes:
    normal           - Serve file normally with Content-Length and Range support
    503_then_ok      - Return 503 for first N requests, then serve normally
    429_then_ok      - Return 429 for first N requests, then serve normally
    drop_mid_transfer - Send first half of file then close connection, N times
    404              - Always return 404
    stall            - Send headers + a few bytes then sleep forever
"""

import argparse
import os
import signal
import sys
import threading
import time
from http.server import HTTPServer, BaseHTTPRequestHandler


class MockHandler(BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.1"  # Use HTTP/1.1 for proper keep-alive
    file_data = b""
    mode = "normal"
    fail_count = 2
    request_count = 0
    lock = threading.Lock()

    def log_message(self, format, *args):
        # Suppress request logging to avoid polluting test output
        pass

    def do_GET(self):
        with MockHandler.lock:
            MockHandler.request_count += 1
            req_num = MockHandler.request_count

        mode = MockHandler.mode
        data = MockHandler.file_data
        fail_count = MockHandler.fail_count

        # Parse Range header
        range_start = 0
        range_header = self.headers.get("Range")
        if range_header and range_header.startswith("bytes="):
            range_spec = range_header[6:]
            if range_spec.endswith("-"):
                range_start = int(range_spec[:-1])
            elif "-" in range_spec:
                parts = range_spec.split("-")
                range_start = int(parts[0])

        if mode == "404":
            self.send_response(404)
            self.send_header("Content-Length", "0")
            self.end_headers()
            return

        if mode == "503_then_ok":
            if req_num <= fail_count:
                self.send_response(503)
                self.send_header("Content-Length", "0")
                self.end_headers()
                return

        if mode == "429_then_ok":
            if req_num <= fail_count:
                self.send_response(429)
                self.send_header("Content-Length", "0")
                self.end_headers()
                return

        if mode == "drop_mid_transfer":
            if req_num <= fail_count:
                # Send headers indicating full length, but only send half
                remaining = data[range_start:]
                half = len(remaining) // 2
                if half < 1:
                    half = 1
                if range_start > 0:
                    self.send_response(206)
                    self.send_header("Content-Range",
                                     "bytes %d-%d/%d" % (range_start,
                                                         len(data) - 1,
                                                         len(data)))
                    self.send_header("Content-Length", str(len(remaining)))
                else:
                    self.send_response(200)
                    self.send_header("Content-Length", str(len(data)))
                self.send_header("Accept-Ranges", "bytes")
                self.send_header("Connection", "close")
                self.end_headers()
                try:
                    self.wfile.write(remaining[:half])
                    self.wfile.flush()
                except BrokenPipeError:
                    pass
                # Force-close the socket to simulate a drop
                import socket
                try:
                    self.connection.setsockopt(
                        socket.SOL_SOCKET, socket.SO_LINGER,
                        bytes([1, 0, 0, 0, 0, 0, 0, 0]))
                    self.connection.shutdown(socket.SHUT_RDWR)
                except Exception:
                    pass
                try:
                    self.connection.close()
                except Exception:
                    pass
                return

        if mode == "stall":
            # Send headers and a few bytes, then sleep forever
            self.send_response(200)
            self.send_header("Content-Length", str(len(data)))
            self.send_header("Accept-Ranges", "bytes")
            self.end_headers()
            try:
                self.wfile.write(data[:min(10, len(data))])
                self.wfile.flush()
                # Sleep until killed
                while True:
                    time.sleep(3600)
            except (BrokenPipeError, ConnectionResetError):
                pass
            return

        # Normal serving (with Range support)
        remaining = data[range_start:]
        if range_start > 0:
            self.send_response(206)
            self.send_header("Content-Range",
                             "bytes %d-%d/%d" % (range_start,
                                                 len(data) - 1,
                                                 len(data)))
            self.send_header("Content-Length", str(len(remaining)))
        else:
            self.send_response(200)
            self.send_header("Content-Length", str(len(data)))
        self.send_header("Accept-Ranges", "bytes")
        self.end_headers()
        self.wfile.write(remaining)


def main():
    parser = argparse.ArgumentParser(description="Mock HTTP server for testing")
    parser.add_argument("--mode", required=True,
                        choices=["normal", "503_then_ok", "429_then_ok",
                                 "drop_mid_transfer", "404", "stall"])
    parser.add_argument("--file", required=True, help="File to serve")
    parser.add_argument("--fail-count", type=int, default=2,
                        help="Number of requests to fail before succeeding")
    parser.add_argument("--port", type=int, default=0,
                        help="Port to listen on (0 for auto)")
    args = parser.parse_args()

    with open(args.file, "rb") as f:
        MockHandler.file_data = f.read()

    MockHandler.mode = args.mode
    MockHandler.fail_count = args.fail_count

    server = HTTPServer(("127.0.0.1", args.port), MockHandler)
    port = server.server_address[1]

    # Print port to stdout so the test program can read it
    print(port, flush=True)

    # Handle SIGTERM gracefully
    def handle_sigterm(signum, frame):
        server.shutdown()
        sys.exit(0)

    signal.signal(signal.SIGTERM, handle_sigterm)

    server.serve_forever()


if __name__ == "__main__":
    main()
