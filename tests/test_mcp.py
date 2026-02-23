import sys
from mcp import ClientSession, StdioServerParameters
from mcp.client.sse import sse_client
import asyncio

async def run():
    if len(sys.argv) < 2:
        print("Usage: python test_mcp.py <port>")
        sys.exit(1)
        
    port = sys.argv[1]
    async with sse_client(url=f"http://localhost:{port}/mcp/sse") as streams:
        async with ClientSession(streams[0], streams[1]) as session:
            await session.initialize()
            # List Prompts
            prompts = await session.list_prompts()
            print("--- MCP Prompts ---")
            for prompt in prompts.prompts:
                print(f"Name: {prompt.name}")
                print(f"Description: {prompt.description}")
                print("-" * 20)
                
                print(f"\n--- Fetching {prompt.name} ---")
                try:
                    p = await session.get_prompt(prompt.name)
                    for msg in p.messages:
                        # Depending on the SDK structure, the content text usually sits under msg.content.text
                        if hasattr(msg.content, 'text'):
                            print(msg.content.text[:200] + "... [TRUNCATED]")
                        else:
                            print(str(msg.content)[:200] + "... [TRUNCATED]")
                except Exception as e:
                    print(f"Error fetching prompt: {e}")

            # List Resources
            resources = await session.list_resources()
            print("\n--- MCP Resources ---")
            for resource in resources.resources:
                print(f"Name: {resource.name}")
                print(f"URI: {resource.uri}")
                print("-" * 20)

            # Read SPARQL Schema Resource
            print("\n--- Reading sparql://schema ---")
            try:
                content = await session.read_resource("sparql://schema")
                print(content.contents[0].text)
            except Exception as e:
                print(f"Error reading resource: {e}")

if __name__ == "__main__":
    asyncio.run(run())
