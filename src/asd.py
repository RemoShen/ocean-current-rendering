import asyncio
import websockets

async def websocket_handler(websocket, path):
    while True:
        message = await websocket.recv()
        await websocket.send(f"Received: {message}")

start_server = websockets.serve(websocket_handler, "localhost", 8080)

asyncio.get_event_loop().run_until_complete(start_server)
asyncio.get_event_loop().run_forever()
