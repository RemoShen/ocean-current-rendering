export function parseNc(path: string) {
    const socket = new WebSocket('ws://localhost:8080');
    socket.addEventListener('open', function (event) {
        console.log('Connected to WS Server');
        socket.send('Hello Server!');
    }
    );
    socket.addEventListener('message', function (event) {
        console.log('Message from server ', event.data);
    }
    );
    socket.addEventListener('close', function (event) {
        console.log('Disconnected from WS Server');
    }
    );
}